import csv
import logging
import multiprocessing
import pickle  # S403
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import Mol, rdSubstructLibrary
from tqdm import tqdm

from chemistry_helpers import process_smiles
from sdf_helpers import find_structures_bytes_ranges, mmap_file, write_mols_to_sdf

RDLogger.DisableLog("rdApp.*")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def export_descriptors_to_csv(descriptors: dict, path: Path) -> None:
    file_exists = path.exists()
    with open(path, "a") as f:
        headers = ["smiles"] + list(descriptors[next(iter(descriptors))].keys())
        csv_writer = csv.writer(f)
        if not file_exists:
            csv_writer.writerow(headers)
        for structure, properties in descriptors.items():
            row = [structure] + [properties.get(key, "") for key in headers[1:]]
            csv_writer.writerow(row)


def load_processed_smiles(path: Path) -> set:
    processed_smiles_set: set = set()
    processed_smiles_file = path / "smiles_processed.csv"
    if processed_smiles_file.exists():
        with open(processed_smiles_file) as f:
            reader = csv.reader(f)
            try:
                next(reader)  # Skip header
            except StopIteration:
                return processed_smiles_set
            for row in reader:
                if len(row) >= 2:
                    processed_smiles_set.add(
                        row[1]
                    )  # cleaned smiles from previous runs
    return processed_smiles_set


def _quick_canonical(smiles: str) -> str | None:
    """Lightweight canonicalization to approximate standardized output for duplicate filtering.

    This is cheaper than full `process_smiles` standardization and helps avoid
    reprocessing between runs when raw input differs but canonical form matches
    a previously processed structure.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        return Chem.MolToSmiles(mol)
    except Exception:
        return None


def run(path: Path) -> None:
    processed_smiles_set = load_processed_smiles(path)
    # Track duplicates prevented within this batch (both raw and quick canonical forms)
    batch_seen: set[str] = set()
    smileses: list[str] = []
    links: list[int] = []
    duplicates_skipped = 0
    canonical_hits_skipped = 0

    structures_csv = path / "structures.csv"
    if not structures_csv.exists():
        logging.error("structures.csv not found at %s", structures_csv)
        return

    with open(structures_csv) as f:
        reader = csv.reader(f)
        try:
            next(reader)
        except StopIteration:
            logging.warning("structures.csv empty; nothing to do")
            return
        for row in reader:
            if len(row) < 3:
                continue
            c, smi, cano = row
            raw = smi.strip() if smi else cano.strip()
            if not raw:
                continue
            # Quick canonical form for filtering; fall back to raw if canonicalization fails
            canon = _quick_canonical(raw) or raw
            # Skip if we already processed this (cleaned) earlier run OR in this batch
            if canon in processed_smiles_set:
                canonical_hits_skipped += 1
                continue
            if canon in batch_seen:
                duplicates_skipped += 1
                continue
            batch_seen.add(canon)
            smileses.append(
                raw
            )  # keep original (will be cleaned later in process_smiles)
            links.append(int(c))

    if not smileses:
        logging.info("No new SMILES to process. Exiting.")
        return

    logging.info(
        "Queued %d new SMILES (skipped %d intra-batch duplicates, %d previously processed)",
        len(smileses),
        duplicates_skipped,
        canonical_hits_skipped,
    )

    max_workers = multiprocessing.cpu_count()

    smis = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
    fps = rdSubstructLibrary.PatternHolder()
    library = rdSubstructLibrary.SubstructLibrary(smis, fps)

    mols_h = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
    fps_h = rdSubstructLibrary.PatternHolder()
    library_h = rdSubstructLibrary.SubstructLibrary(mols_h, fps_h)

    inchis: list[str] = []
    inchikeys: list[str] = []
    formulas: list[str] = []
    sdf_blocks: list[tuple[int, str]] = []
    p_smileses: list[str] = []
    p_smols = []
    p_sim_fps = []
    p_sim_h_fps = []
    p_links: list[int] = []
    smis_no_stereo: list[str] = []
    inchis_no_stereo: list[str] = []
    inchikeys_no_stereo: list[str] = []
    descriptors_r: dict = {}

    logging.info("Generating the chemical libraries (workers=%d)", max_workers)
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Enumerate to retain mapping index -> structure ID (links list)
        futures = {
            executor.submit(process_smiles, item): item for item in enumerate(smileses)
        }
        for future in tqdm(
            as_completed(futures),
            total=len(smileses),
            desc="Processing SMILES",
        ):
            try:
                result = future.result()
            except Exception:
                logging.exception("Worker future failed")
                continue
            if not result:
                continue
            (
                nid,
                smiles_original,
                smol,
                smiles_clean,
                inchi_clean,
                inchikey_clean,
                formulas_clean,
                mol_block,
                sim_fp,
                sub_fp,
                desc_rdkit,
                mol_h,
                sim_fp_h,
                sub_fp_h,
                smiles_no_stereo_val,
                inchi_no_stereo_val,
                inchikey_no_stereo_val,
            ) = result

            # Hydrogens variant
            mols_h.AddMol(Mol(mol_h))
            fps_h.AddFingerprint(sub_fp_h)
            p_sim_h_fps.append(sim_fp_h)

            # Base variant
            smis.AddSmiles(smiles_clean)
            fps.AddFingerprint(sub_fp)
            p_sim_fps.append(sim_fp)

            inchis.append(inchi_clean)
            inchikeys.append(inchikey_clean)
            formulas.append(formulas_clean)
            sdf_blocks.append((links[nid], mol_block))
            p_smols.append(smol)
            # TODO check which one to use
            # See #109
            # p_smileses.append(smiles)
            p_smileses.append(smiles_clean)
            smis_no_stereo.append(smiles_no_stereo_val)
            inchis_no_stereo.append(inchi_no_stereo_val)
            inchikeys_no_stereo.append(inchikey_no_stereo_val)
            descriptors_r[smiles_clean] = desc_rdkit
            p_links.append(links[nid])

    logging.info("Finished generating the chemical libraries")

    logging.info("Generating and exporting SDF")
    write_mols_to_sdf(path / "lotus.sdf", sdf_blocks)

    logging.info("Indexing SDF")
    mapped_sdf = mmap_file(path / "lotus.sdf")
    structures_ranges = find_structures_bytes_ranges(mapped_sdf)

    logging.info("Generating database")
    database = {
        "structure_wid": p_links,
        "structure_sim_fps": p_sim_fps,
        "structure_sim_h_fps": p_sim_h_fps,
        "structure_library": library.Serialize(),
        "structure_library_h": library_h.Serialize(),
        "structure_id": {i[1]: i[0] for i in enumerate(p_links)},
    }

    logging.info("Exporting primary table")
    smiles_file_path = path / "structures_table.csv"
    file_exists = smiles_file_path.exists()
    with open(smiles_file_path, "a") as f:  # Append mode
        csv_file = csv.writer(f)
        if not file_exists:
            # TODO check if we want the clean SMILES or not in this table
            csv_file.writerow([
                "structure",
                "structure_smiles",
                "structure_smiles_no_stereo",
                "structure_inchi",
                "structure_inchi_no_stereo",
                "structure_inchikey",
                "structure_inchikey_no_stereo",
                "structure_formula",
            ])
        csv_file.writerows(
            zip(
                p_links,
                p_smileses,
                smis_no_stereo,
                inchis,
                inchis_no_stereo,
                inchikeys,
                inchikeys_no_stereo,
                formulas,
                strict=False,
            ),
        )

    # TODO this is not finished.
    logging.info("Exporting blocks table")
    blocks_path = path / "structures_blocks_table.csv"
    file_exists = blocks_path.exists()
    with open(blocks_path, "a") as f:
        csv_file = csv.writer(f)
        if not file_exists:
            csv_file.writerow(["structure", "block_range"])
        csv_file.writerows(
            zip(structures_ranges, structures_ranges.values(), strict=False)
        )

    logging.info("Exporting rdkit descriptors")
    export_descriptors_to_csv(descriptors_r, path / "descriptors_rdkit.csv")

    logging.info("Exporting processed smiles")
    smiles_processed_path = path / "smiles_processed.csv"
    file_exists = smiles_processed_path.exists()
    with open(smiles_processed_path, "a") as f:
        csv_file = csv.writer(f)
        if not file_exists:
            csv_file.writerow(["structure", "structure_smiles"])
        csv_file.writerows(zip(p_links, p_smileses, strict=False))

    logging.info("Exporting database")
    with open(path / "database_chemo.pkl", "wb") as f:
        pickle.dump(database, f)

    logging.info("Finished exporting")


if __name__ == "__main__":
    run(Path("data"))

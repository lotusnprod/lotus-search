import mmap
from collections import deque
from collections.abc import Iterable
from pathlib import Path

from rdkit.Chem import MolFromMolBlock, SDWriter


def mmap_file(path: Path) -> mmap.mmap:
    with open(path) as file:
        mmapped_file = mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ)
    return mmapped_file


def find_structures_bytes_ranges(mmapped_file: mmap.mmap) -> dict[int, tuple[int, int]]:
    structures_ranges: dict[int, tuple[int, int]] = {}
    start_offset: int = 0
    prev_line: bytes = b""
    prev_prev_line: bytes = b""
    for line in iter(mmapped_file.readline, b""):
        if line.startswith(b"$$$$"):
            end_offset = mmapped_file.tell() - len(line)
            identifier = int(prev_prev_line.strip().decode())
            structures_ranges[identifier] = (start_offset, end_offset)
            start_offset = mmapped_file.tell()

        prev_prev_line, prev_line = prev_line, line
    return structures_ranges


def read_selected_ranges(mmapped_file: mmap.mmap, ranges_to_read: list[tuple[int, int]]) -> str:
    selected_lines: deque[str] = deque()

    for start, end in ranges_to_read:
        selected_lines.append(mmapped_file[start:end].decode())

    return "".join(selected_lines)


def write_mols_to_sdf(path: Path, sdf_blocks: Iterable[tuple[int, str]]) -> None:
    with open(str(path), "a") as f:
        sorted_sdf_blocks = sorted(sdf_blocks, key=lambda x: x[0])
        with SDWriter(f) as w:
            for wid, sdf_block in sorted_sdf_blocks:
                mol = MolFromMolBlock(sdf_block)
                if mol:
                    mol.SetProp("WID", str(wid))
                    w.write(mol)


def main():
    file_path = "data/lotus.sdf"
    mmapped_file = mmap_file(file_path)
    structures_ranges = find_structures_bytes_ranges(mmapped_file)
    # print(structures_ranges)
    start_time = time.time()
    ranges_to_read = [structures_ranges[key] for key in list(structures_ranges.keys())[:100000]]
    if ranges_to_read:
        selected_lines = read_selected_ranges(mmapped_file, ranges_to_read)
    else:
        logging.info("No '$$$$' occurrences found in the file.")
    end_time = time.time()
    logging.info(f"Time taken to get the blocks: {end_time - start_time} seconds")
    # print(selected_lines)


if __name__ == "__main__":
    import logging
    import time

    logging.basicConfig(level=logging.INFO)
    main()

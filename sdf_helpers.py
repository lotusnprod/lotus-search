import mmap

from pathlib import Path
from collections import deque
from rdkit.Chem import MolFromMolBlock, SDWriter


def mmap_file(file_path):
    with open(file_path, "r") as file:
        mmapped_file = mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ)
    return mmapped_file


def find_structures_line_ranges(mmapped_file):
    structures_ranges = {}
    start_line = 1
    prev_line = prev_prev_line = b""
    for line_number, line in enumerate(iter(mmapped_file.readline, b""), start=1):
        if line.startswith(b"$$$$"):
            end_line = line_number
            identifier = int(prev_prev_line.strip().decode())
            structures_ranges[identifier] = (start_line, end_line)
            start_line = line_number + 1
        prev_prev_line, prev_line = prev_line, line
    return structures_ranges


def read_selected_ranges(mmapped_file, ranges_to_read):
    selected_lines = deque()
    for start, end in ranges_to_read:
        selected_lines.append(mmapped_file[start:end].decode())
    return selected_lines


def write_mols_to_sdf(path: Path, sdf_blocks):
    with SDWriter(str(path)) as w:
        for wid, sdf_block in sdf_blocks:
            mol = MolFromMolBlock(sdf_block)
            if mol:
                mol.SetProp("WID", str(wid))
                w.write(mol)


def main():
    file_path = "data/lotus.sdf"
    mmapped_file = mmap_file(file_path)
    structures_ranges = find_structures_line_ranges(mmapped_file)
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
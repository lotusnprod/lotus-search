import time
import logging
import mmap
from collections import deque


def find_structures_line_ranges(file_path):
    structures_ranges = {}
    with open(file_path, "r") as file:
        mmapped_file = mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ)
        start_line = 1
        prev_line = prev_prev_line = b''
        for line_number, line in enumerate(iter(mmapped_file.readline, b""), start=1):
            if line.startswith(b"$$$$"):
                end_line = line_number
                identifier = prev_prev_line.strip().decode()  # Use the value two lines above as the identifier
                structures_ranges[identifier] = (start_line, end_line)
                start_line = line_number + 1
            prev_prev_line, prev_line = prev_line, line
    return structures_ranges, mmapped_file


def read_selected_ranges(mmapped_file, ranges_to_read):
    selected_lines = deque()
    for start, end in ranges_to_read:
        selected_lines.append(mmapped_file[start:end].decode())
    return selected_lines


def main():
    file_path = "data/lotus.sdf"
    structures_ranges, mmapped_file = find_structures_line_ranges(file_path)
    # print(structures_ranges)

    # Admitting the indices are already available
    start_time = time.time()
    # TODO this will allow to get the ranges from a (list of) SID(s)
    ranges_to_read = [structures_ranges[key] for key in list(structures_ranges.keys())[:100000]]

    if ranges_to_read:
        selected_lines = read_selected_ranges(mmapped_file, ranges_to_read)
    else:
        logging.info("No '$$$$' occurrences found in the file.")

    end_time = time.time()
    logging.info(f"Time taken to get the blocks: {end_time - start_time} seconds")
    # print(selected_lines)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
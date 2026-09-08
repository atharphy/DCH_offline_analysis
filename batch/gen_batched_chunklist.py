import sys
from collections import defaultdict

BATCH_SIZE = 15

def main():
    groups = defaultdict(list)
    with open("batch/chunklist.txt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            year, idx, chunk = line.split(",")
            groups[(year, int(idx))].append(int(chunk))

    with open("batch/batched_chunklist.txt", "w") as out:
        for (year, idx), chunks in groups.items():
            chunks.sort()
            for start in range(0, len(chunks), BATCH_SIZE):
                batch = chunks[start:start + BATCH_SIZE]
                out.write(f"{year},{idx},{batch[0]},{len(batch)}\n")

if __name__ == "__main__":
    main()

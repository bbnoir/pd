#!/usr/bin/env python3

import subprocess
import os
import random
from concurrent.futures import ThreadPoolExecutor, as_completed

def main():
    """
    Main function to run all test targets.
    """
    # Ensure the script is run from the correct directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    targets = ["0", "1", "2", "3", "4", "5"]

    target_scores = []
    for target in targets:
        make_target = f"t{int(target)+1}"
        result = subprocess.run(
            ["make", make_target],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        target_score = float(result.stdout.split()[-1])
        target_scores.append(target_score)
        print(result.stdout)

    avg_score = sum(target_scores) / len(target_scores)
    print(f"Average score: {avg_score:.3f}")
    print()

if __name__ == "__main__":
    main()
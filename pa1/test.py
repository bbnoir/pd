#!/usr/bin/env python3

import subprocess
import os
import random
from concurrent.futures import ThreadPoolExecutor, as_completed

def run_with_seed(target, seed):
    try:
        # print(f"Running test for input_{target} with seed {seed}")
        input_file = f"input_pa1/input_{target}.dat"
        output_file = f"out/input_{target}.dat"
        result = subprocess.run(
            ["./bin/fm", input_file, output_file, str(seed)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        # print(result.stdout)
        if result.stderr:
            print(f"Warnings/Errors:\n{result.stderr}")
        else:
            return int(result.stdout.strip())
    except subprocess.CalledProcessError as e:
        print(f"Error while running target {target}:\n{e.stderr}")

def main():
    """
    Main function to run all test targets.
    """
    # Ensure the script is run from the correct directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    targets = ["0", "1", "2", "3", "4", "5"]
    # targets = ["1", "2"]

    gen_seed = 373
    random.seed(gen_seed)
    random_num = 1000
    seeds = [random.randint(0, 100000) for _ in range(random_num)]

    best5_seeds = []
    best5_scores = []

    # Run each test target
    for target in targets:
        print(f"Running test for input_{target}")
        scores_and_seeds = []
        
        with ThreadPoolExecutor(max_workers=16) as executor:
            future_to_seed = {executor.submit(run_with_seed, target, seed): seed for seed in seeds}
            
            for future in as_completed(future_to_seed):
                seed = future_to_seed[future]
                result = future.result()
                if result is not None:
                    # print(f"seed {seed} -> {result}")
                    scores_and_seeds.append((result, seed))
                else:
                    print(f"[ERROR] Failed to run test for input_{target} with seed {seed}")
        
        # Sort by score (ascending) and get top 5
        scores_and_seeds.sort()
        top_5 = scores_and_seeds[:5]
        
        for rank, (score, seed) in enumerate(top_5, 1):
            make_target = f"t{int(target)+1}"
            seed_arg = f"SEED={seed}"
            result = subprocess.run(
                ["make", make_target, seed_arg],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            # print(result.stdout)
            final_score = result.stdout.split()[-1]
            print(f"Seed {seed} -> Cut {score} -> Final score {final_score}")
            best5_scores.append(float(final_score))
            best5_seeds.append(seed)

if __name__ == "__main__":
    main()
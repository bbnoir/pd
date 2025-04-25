#!/usr/bin/env python3

import subprocess
import os
import random
from concurrent.futures import ThreadPoolExecutor, as_completed

input_dir = "input_pa2"
# targets = ["ami33"]
targets = ["ami33", "ami49", "apte", "hp", "xerox"]
test_dir = "testdir"
# subdir
output_dir = f"{test_dir}/out"
log_dir = f"{test_dir}/log"

executable = "./bin/fp"
evaluator = "./evaluator/evaluator.sh"
alpha = 0.5

def run_with_seed(target, seed):
    try:
        # print(f"Running test for input_{target} with seed {seed}")
        block_file = f"{input_dir}/{target}.block"
        nets_file = f"{input_dir}/{target}.nets"
        output_file = f"{output_dir}/{target}_{seed}.out"
        command = [executable, str(alpha), block_file, nets_file, output_file, str(seed)]
        result = subprocess.run(
            command,
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

    if not os.path.exists(test_dir):
        os.makedirs(test_dir)
        os.makedirs(output_dir)
        os.makedirs(log_dir)

    gen_seed = 8383
    random.seed(gen_seed)
    random_num = 2000
    seeds = [random.randint(0, 100000) for _ in range(random_num)]

    highest_score = []

    # Run each test target
    for target in targets:
        print(f"Running test for {target}")
        log_file = f"{log_dir}/{target}.log"
        with open(log_file, "w") as log:
            log.write(f"Running test for {target}\n")
            scores_and_seeds = []

            with ThreadPoolExecutor(max_workers=16) as executor:
                future_to_seed = {executor.submit(run_with_seed, target, seed): seed for seed in seeds}
                
                for future in as_completed(future_to_seed):
                    seed = future_to_seed[future]
                    result = future.result()
                    if result is not None:
                        # log.write(f"seed {seed} -> {result}\n")
                        if result != -1:
                            scores_and_seeds.append((result, seed))
                    else:
                        log.write(f"[ERROR] Failed to run test for input_{target} with seed {seed}\n")
            
            # Sort by score (ascending) and get top 5
            scores_and_seeds.sort()
            top_5 = scores_and_seeds[:10]

            highest = -100.0
            
            log.write(f"Top 10 scores for {target}:\n")
            for rank, (score, seed) in enumerate(top_5, 1):
                block_file = f"{input_dir}/{target}.block"
                nets_file = f"{input_dir}/{target}.nets"
                out_file = f"{output_dir}/{target}_{seed}.out"
                command = ["bash", evaluator, block_file, nets_file, out_file, str(alpha)]
                result = subprocess.run(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True
                )
                # log.write(result.stdout)
                final_score = result.stdout.split()[-1]
                highest = max(highest, float(final_score))
                log.write(f"Seed {seed} -> Cost {score} -> Final score {final_score}\n")
            
            highest_score.append((target, highest))
        
        # clean out
        for seed in seeds:
            out_file = f"{output_dir}/{target}_{seed}.out"
            if os.path.exists(out_file):
                os.remove(out_file)
        
    print("Highest scores:")
    for target, score in highest_score:
        print(f"{target}: {score}")
    print(f"Average score: {sum(score for _, score in highest_score) / len(highest_score):.3f}")
    

if __name__ == "__main__":
    main()
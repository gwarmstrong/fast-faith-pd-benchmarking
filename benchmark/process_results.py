import os
import json
from util import safe_dir

def process_dir(input_dir, output_file=None):
    # TODO: sanity check on directory
    file_names = [os.path.join(input_dir, file_) for file_ in os.listdir(input_dir)]

    all_results = [process_results_file(file_) for file_ in file_names]

    if output_file is None:
        output_path = 'results.json'
    else:
        #safe_dir(output_file)
        output_path = os.path.join(output_file)
    
    with open(output_path, 'w') as fp:
        json.dump(all_results, fp)
    

def process_results_file(results_file):
    with open(results_file, 'r') as fp:
        lines = fp.readlines()

    lines = [line.strip() for line in lines]

    results = {}

    for i in range(4):
        this_line = lines.pop(0)
        k, v = this_line.split('-', 1)
        results[k] = v

    results['method'] = lines[0].split('-')[1]
    results['python_time'] = lines[1].split('-')[1]

    times = filter(lambda x: len(x) > 0, lines[1].split('  '))
    for time in times:
        t, k = time.split(' ')
        results[k] = t

    for line in lines[3:]:
        t, k = line.split('  ')
        results[k] = t

    return results


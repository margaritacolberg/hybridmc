import asyncio
import csv


async def find_paths(start_bitstring, termination_bitstring, current_path, all_paths):
    if start_bitstring == termination_bitstring:
        all_paths.append(current_path.copy())
        return

    for i in range(len(start_bitstring)):
        if start_bitstring[i] != termination_bitstring[i]:
            next_bitstring = start_bitstring[:i] + \
                             str(int(not int(start_bitstring[i]))) + \
                             start_bitstring[i + 1:]
            current_path.append(next_bitstring)
            await find_paths(next_bitstring, termination_bitstring, current_path, all_paths)
            current_path.pop()


async def process_path(path, transition_info):
    summed_entropy = 0
    summed_weights = 0
    for i in range(len(path) - 1):
        bitstring_in = path[i]
        bitstring_out = path[i + 1]
        for transition in transition_info:
            if transition[0] == bitstring_in and transition[1] == bitstring_out:
                summed_entropy += transition[2]
                summed_weights += transition[3] ** 2
                break
    return summed_entropy, summed_weights


async def process_all_paths(all_paths, transition_info):
    master_entropy_weights_list = []
    tasks = []
    for path in all_paths:
        task = asyncio.create_task(process_path(path, transition_info))
        tasks.append(task)
    results = await asyncio.gather(*tasks)
    for result in results:
        master_entropy_weights_list.append(result)
    return master_entropy_weights_list


async def main():
    start_bitstring = "00"
    termination_bitstring = "11"
    all_paths = []
    await find_paths(start_bitstring, termination_bitstring, [start_bitstring], all_paths)

    transition_info = []  # Load transition info from CSV file

    # Simulated transition info, replace this with actual loading from CSV
    with open('diff_s_bias_with_error.csv', newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            transition_info.append((row[0], row[1], float(row[2]), float(row[3])))

    master_entropy_weights_list = await process_all_paths(all_paths, transition_info)

    print("Master list of (summed_entropy, summed_weights) tuples for each path:")
    for item in master_entropy_weights_list:
        print(item)


if __name__ == "__main__":
    asyncio.run(main())

with open('data.txt', 'r') as source, open('d.txt', 'w') as target:
    for line in source:
        if len(line) == 0:
            continue
        arr = line.split()
        for index in (0, 1, 2, 3, 4, 5, 8):
            target.write(arr[index] + ' ')
        target.write('\n')
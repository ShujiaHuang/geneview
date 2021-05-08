"""
Test some sample straitage
"""
from __future__ import division

def sample1():
    """
    real    0m0.397s
    user    0m0.321s
    sys 0m0.059s
    """
    row_num = 0
    with open('test.fastq', 'r') as I, open('output.fastq', 'w') as O:
        for row1 in I:
            row2 = I.next()
            row3 = I.next()
            row4 = I.next()
            if row_num % 10 == 0:
                O.write(row1)
                O.write(row2)
                O.write(row3)
                O.write(row4)
            row_num += 1

def sample2():
    """
    real    0m0.800s
    user    0m0.682s
    sys 0m0.096s
    """
    import numpy as np

    with open('test.fastq', 'r') as I, open('output.fastq', 'w') as O:
        for row1 in I:
            row2 = I.next()
            row3 = I.next()
            row4 = I.next()
            if np.random.randint(10) == 0:
                O.write(row1)
                O.write(row2)
                O.write(row3)
                O.write(row4)


def sample3():
    """
    """
    import numpy as np

    percent = 30
    with open('test.fastq', 'r') as I, open('output.fastq', 'w') as O:
        for row1 in I:
            row2 = I.next()
            row3 = I.next()
            row4 = I.next()
            if np.random.randint(1, 101) <= percent:
                O.write(row1)
                O.write(row2)
                O.write(row3)
                O.write(row4)

def sample4():
    #from __future__ import division
    import numpy as np

    num_to_sample = 30000
    with open('test.fastq', 'r') as I:
        line_num = sum([1 for r in I])
    total_records = line_num // 4

    percent = (num_to_sample / total_records) * 100
    print line_num, total_records, num_to_sample, percent
    with open('test.fastq', 'r') as I, open('output.fastq', 'w') as O:
        for row1 in I:
            row2 = I.next()
            row3 = I.next()
            row4 = I.next()
            if np.random.randint(1, 101) <= percent:
                O.write(row1)
                O.write(row2)
                O.write(row3)
                O.write(row4)

def sample5():
    import numpy as np

    num_to_sample = 30000
    with open('test.fastq', 'r') as I:
        line_num = sum([1 for r in I])

    total_records = range(line_num // 4 + 1)
    np.random.shuffle(total_records)
    record2keep = set(total_records[:num_to_sample])

    print record2keep
    record_num = 0
    with open('test.fastq', 'r') as I, open('output.fastq', 'w') as O:
        for row1 in I:
            row2 = I.next()
            row3 = I.next()
            row4 = I.next()
            if record_num in record2keep:
                O.write(row1)
                O.write(row2)
                O.write(row3)
                O.write(row4)
            record_num += 1

if __name__ == '__main__':
    sample5()

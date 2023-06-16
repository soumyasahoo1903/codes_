from sklearn.linear_model import ElasticNet
import tsutil as tsutil
import numpy as np
import multiprocessing as mp


mb_np, label_np, mb_TC_np, label_TC_np, mb_AB_np, label_AB_np = tsutil.read_ts_data(1.5)


def calculate_score(mb_np_a, label_np_a, i, j):
    y = (mb_np_a[i, :]/np.linalg.norm(mb_np_a[i, :])).reshape(-1, 1)
    x = label_np_a[j, :]/np.linalg.norm(label_np_a[j, :])
    x = x*1000
    y = y*1000
    regressor = ElasticNet()
    regressor.fit(y, x)
    return [str([i, j]), regressor.score(y, x)]

def main(a):
    scores = []
    if a == 'o':
        label_np_a = label_np
        mb_np_a = mb_np
    if a == 'tc':
        label_np_a = label_TC_np
        mb_np_a = mb_TC_np
    if a == 'ab':
        label_np_a = label_AB_np
        mb_np_a = mb_AB_np
    num_cores = mp.cpu_count()  # Get the number of CPU cores
    pool = mp.Pool(processes=num_cores)  # Create a process pool

    results = []
    for i in range(mb_np_a.shape[0]):
        for j in range(label_np_a.shape[0]):
            results.append(pool.apply_async(calculate_score, (mb_np_a, label_np_a, i, j)))

    for result in results:
        tsutil.save_score(str(result.get()), a, "en")

    pool.close()
    pool.join()


if __name__ == '__main__':
    main('tc')
    main('ab')
    main('o')


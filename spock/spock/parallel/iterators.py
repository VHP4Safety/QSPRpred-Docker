from pebble import ProcessPool
from concurrent.futures import TimeoutError

def batched_iter(iterator, batch_size):
    batch = []
    for item in iterator:
        batch.append(item)
        if len(batch) == batch_size:
            yield batch
            batch = []
    yield batch

def parallel_iter(iterator, process_func, n_cpus, *args, timeout=None,   **kwargs):
    with ProcessPool(max_workers=n_cpus) as pool:
        queue = []  # a queue for our current worker async results, a deque would be faster
        done = False
        while not done or queue:
            try:
                # add our next slice to the pool
                queue.append(pool.schedule(process_func, args=[next(iterator)] + list(args), kwargs=kwargs, timeout=timeout))
            except (StopIteration, TypeError) as exp:  # no more data, clear out the slice iterator
                done = True
            # wait for a free worker or until all remaining workers finish
            while queue and (len(queue)  >= n_cpus or done):
                process = queue.pop(0)  # grab a process response from the top

                if process.done() and type(process._exception ) == TimeoutError:  # check if job timed out
                    yield  process._exception
                    break  # make sure to pop the next item in the queue

                try:
                    # check if result available
                    result = process.result(timeout=0.01)
                    yield result  # yield the result if available
                except TimeoutError:
                    queue.append(process)  # add it back to the queue
                except Exception as exp:
                    yield process._exception
                    break  # make sure to pop the next item in the queue

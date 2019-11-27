import numpy as np
import time

#def foo():
#    print("foo")
#
#def bar(func):
#    func()
#
#bar(foo)

def display_time(func):   #decorator
    def wrapper(*args,**kwargs):
        t_start = time.time()
        result = func(*args,**kwargs)
        t_end = time.time()
        print("Total time: {:.4} s".format(t_end - t_start))
        return result
    return wrapper

@display_time    # call decorator
def sum_n(minnum,maxnum):
    sum_ = 0
    count = 0
    for i in range(minnum,maxnum):
        sum_ = sum_ + i
        if i > 0:
            count = count + 1
    return count

    
count = sum_n(-10,10)
print(count)


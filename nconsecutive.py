##two functions that will find the longest set of consecutive integers in a list/array

def sequ(s):
    rl = {}
    best_range = xrange(0)
    for x in s:
        run = rl[x] = rl.get(x-1, 0) + 1
        r = xrange(x-run+1, x+1)
        if len(r) > len(best_range):
            best_range = r
    return rl, best_range

def sequ2(s):
    rl={}
    maxrun=-1
    for x in s:
        run = rl[x] = rl.get(x-1, 0) + 1
       # print x-run+1, 'to', x
        if run > maxrun:
            maxend, maxrun = x, run
    return rl, range(maxend-maxrun+1, maxend+1)

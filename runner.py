import sys
from main import *

def run(nodes, q, iters, d, times):
    for i in range(int(times)):
        main(nodes, q, iters, d)

if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
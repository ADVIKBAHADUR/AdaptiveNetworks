import sys
from main import *

avg_errors = []
qs = []
avg_c = 0
avg_cs = []

def run(nodes, iters, d, times):
    ratio = 0.5
    q = 0
    for i in range(int(times)):
        avg_error, avg_c = main(nodes, q, iters, d, ratio)
        qs.append(q)
        avg_cs.append(avg_c)
        avg_errors.append(avg_error)
        q = q + 0.02
    fig, ax = plt.subplots()
    print(avg_errors)
    print(avg_cs)
    color = 'tab:blue'
    ax.plot(qs, avg_errors)
    ax.set_ylabel("avg_errors")
    ax.set_xlabel("q")
    plt.show()
    # for i in range(int(times)):
    #     ratios.append(ratio)
    #     avg_error, avg_c_degree = main(nodes, q, iters, d, ratio)
    #     avg_errors.append(avg_error)
    #     ratio = ratio + 0.01
    # print(avg_errors)
    # fig, ax = plt.subplots()
    # color = 'tab:blue'
    # ax.plot(ratios, avg_errors)
    # ax.set_ylabel("Average Errors")
    # ax.set_xlabel('Ratio Central/Outer')
    # ax.tick_params(axis='y', labelcolor=color)

    # ax2 = ax.twinx()

    # color = 'tab:red'
    # ax2.set_ylabel('Average degree of control layer', color=color)
    # ax2.plot(ratio, avg_c_degree, color=color)
    # ax2.tick_params(axis='y', labelcolor=color)
    # fig.tight_layout()

    # plt.show()

if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
import matplotlib.pyplot as plt
import csv
import pdb
from typing import Optional


def read_csv_data(filename: str) -> dict:
    index2key = {}
    result = {}
    with open(filename, encoding="utf8") as f:
        for i, row in enumerate(csv.reader(f)):
            if i == 0:
                for col_num, entry in enumerate(row):
                    index2key[col_num] = entry
                    result[entry] = []
            else:
                for col_num, entry in enumerate(row):
                    result[index2key[col_num]] += [float(entry)]
    return result


def plot_data(title: str, legends: list[str], x_label: str, y_label, data_x: list, data_y: list, boundary: Optional[tuple] = None):
    assert len(legends) == len(data_x) and len(data_x) == len(data_y)
    _fig, ax = plt.subplots()
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    for label, x, y in zip(legends, data_x, data_y):
        ax.plot(x, y, label=label)
    ax.legend()
    if boundary is not None:
        ax.set_ybound(boundary[0], boundary[1])
    plt.savefig("plot.pdf")
    plt.show()


def plot_same_label_csvs(title: str, filenames: list[str], legends: list[str], x_label, y_label, boundary: Optional[tuple] = None):
    data_x = []
    data_y = []
    for file in filenames:
        data = read_csv_data(file)
        data_x += [data[x_label]]
        data_y += [data[y_label]]
    plot_data(title, legends, x_label, y_label, data_x, data_y, boundary)


def milestone4():
    plot_same_label_csvs("Comparison of different timesteps",
                         ["../data/04/timestep 0.001/data.csv", "../data/04/timestep 0.01/data.csv", "../data/04/timestep 0.04/data.csv"],
                         ["timestep 0.001", "timestep 0.01", "timestep 0.04"],
                         "simulated_time",
                         "total_energy")


def milestone5_6():
    time_5 = read_csv_data("../data/05/time.csv")
    time_6 = read_csv_data("../data/06/time.csv")
    plot_data("Execution time with and without cutoff",
              ["without cutoff", "with cutoff"],
              "number of atoms",
              "execution time [s]",
              [[d ** 3 for d in time_5["dim_x"]], [d ** 3 for d in time_6["dim_x"]]],
              [time_5["time"], time_6["time"]])


def calculate_heat_shit(filename: str, tmp=800.0) -> dict:
    data = read_csv_data(filename)
    index = 0
    for i, t in enumerate(data["temperature"]):
        if t >= tmp:
            index = i
            break
    else:
        raise ValueError("should be able to reach at least 800 Kelvin")
    delta_q = data["total_energy"][index] - data["total_energy"][0]
    delta_tmp = data["temperature"][index] - data["temperature"][0]
    heat_capacity = delta_q / delta_tmp
    # in the energy-tmp plot the curve is linear outside the melting and evaporation point
    total_heat = data["temperature"][-1] - data["temperature"][0]
    total_energy = data["total_energy"][-1] - data["total_energy"][0]
    latent_heat = total_energy - heat_capacity * total_heat
    melting_point = -1.0
    start_energy = data["total_energy"][0]
    start_tmp = data["temperature"][0]
    for i, (tot, tp) in enumerate(zip(data["total_energy"], data["temperature"])):
        d_q = tot - start_energy
        tp_thres = d_q * (1 / heat_capacity) + start_tmp - (latent_heat / 2) * (1 / heat_capacity)
        if tp <= tp_thres:
            melting_point = tp
            break
    return {
        "heat_capacity": heat_capacity,
        "latent_heat": latent_heat,
        "melting_point": melting_point
    }


def milestone7_heat():
    heat_cluster_923 = calculate_heat_shit("../data/07/cluster_923/data.csv")
    heat_cluster_1415 = calculate_heat_shit("../data/07/cluster_1415/data.csv")
    heat_cluster_2057 = calculate_heat_shit("../data/07/cluster_2057/data.csv")
    heat_cluster_2869 = calculate_heat_shit("../data/07/cluster_2869/data.csv")
    heat_cluster_3871 = calculate_heat_shit("../data/07/cluster_3871/data.csv")
    plot_data("Heat capacity of different cluster size",
              ["heat capacity"],
              "cluster size",
              "heat capacity",
              [[923, 1415, 2057, 2869, 3871]],
              [[heat_cluster_923["heat_capacity"], heat_cluster_1415["heat_capacity"],
                heat_cluster_2057["heat_capacity"], heat_cluster_2869["heat_capacity"],
                heat_cluster_3871["heat_capacity"]]])
    plot_data("Latent of different cluster size",
              ["latent heat"],
              "cluster size",
              "latent heat",
              [[923, 1415, 2057, 2869, 3871]],
              [[heat_cluster_923["latent_heat"], heat_cluster_1415["latent_heat"],
                heat_cluster_2057["latent_heat"], heat_cluster_2869["latent_heat"],
                heat_cluster_3871["latent_heat"]]])
    plot_data("Melting point of different cluster size",
              ["melting point"],
              "cluster size",
              "melting point",
              [[923, 1415, 2057, 2869, 3871]],
              [[heat_cluster_923["melting_point"], heat_cluster_1415["melting_point"],
                heat_cluster_2057["melting_point"], heat_cluster_2869["melting_point"],
                heat_cluster_3871["melting_point"]]])


def milestone7_energy():
    cluster_923 = read_csv_data("../data/07/cluster_923/data.csv")
    cluster_1415 = read_csv_data("../data/07/cluster_1415/data.csv")
    cluster_2057 = read_csv_data("../data/07/cluster_2057/data.csv")
    cluster_2869 = read_csv_data("../data/07/cluster_2869/data.csv")
    cluster_3871 = read_csv_data("../data/07/cluster_3871/data.csv")
    cluster_923_x = [x - cluster_923["total_energy"][0] for x in cluster_923["total_energy"]]
    cluster_1415_x = [x - cluster_1415["total_energy"][0] for x in cluster_1415["total_energy"]]
    cluster_2057_x = [x - cluster_2057["total_energy"][0] for x in cluster_2057["total_energy"]]
    cluster_2869_x = [x - cluster_2869["total_energy"][0] for x in cluster_2869["total_energy"]]
    cluster_3871_x = [x - cluster_3871["total_energy"][0] for x in cluster_3871["total_energy"]]
    plot_data("Heating of different cluster sizes",
              ["cluster_923 delta_q 1.5", "cluster_1415 delta_q 2.2", "cluster_2057_x delta_q 3.3", "cluster_2869_x delta_q 4.6", "cluster_3871 delta_q 6"],
              "added energy",
              "temperature",
              [cluster_923_x, cluster_1415_x, cluster_2057_x, cluster_2869_x, cluster_3871_x],
              [cluster_923["temperature"], cluster_1415["temperature"], cluster_2057["temperature"], cluster_2869["temperature"], cluster_3871["temperature"]])


def milestone7():
    milestone7_energy()
    milestone7_heat()


def milestone8():
    plot_same_label_csvs("Energy conservation when working with multiple MPI processes",
                         ["../data/08/1/data.csv", "../data/08/2/data.csv", "../data/08/4/data.csv", "../data/08/8/data.csv"],
                         ["1 process", "2 processes", "4 processes", "8 processes"],
                         "simulated_time",
                         "total_energy",
                         boundary=(-3342, -3344))


def milestone9():
    plot_same_label_csvs("Strain-Stress plot for different temperatures and strain rates (small whisker)",
                         ["../data/09/datasmallt0l0.05.csv",
                          "../data/09/datasmallt0l0.1.csv",
                          "../data/09/datasmallt100l0.05.csv",
                          "../data/09/datasmallt100l0.1.csv",],
                         ["length_increase=0.05 temperature=0", "length_increase=0.1 temperature=0", "length_increase=0.05 temperature=100", "length_increase=0.1 temperature=100",],
                         "strain",
                         "stress")
    plot_same_label_csvs("Strain-Stress plot for different temperatures and strain rates (medium whisker)",
                         ["../data/09/datamediumt0l0.05.csv",
                          "../data/09/datamediumt0l0.1.csv",
                          "../data/09/datamediumt100l0.05.csv",
                          "../data/09/datamediumt100l0.1.csv",],
                         ["length_increase=0.05 temperature=0", "length_increase=0.1 temperature=0", "length_increase=0.05 temperature=100", "length_increase=0.1 temperature=100",],
                         "strain",
                         "stress")
    plot_same_label_csvs("Strain-Stress plot for different temperatures and strain rates (large whisker)",
                         ["../data/09/datalarget0l0.05.csv",
                          "../data/09/datalarget0l0.1.csv",
                          "../data/09/datalarget100l0.05.csv",
                          "../data/09/datalarget100l0.1.csv",],
                         ["length_increase=0.05 temperature=0", "length_increase=0.1 temperature=0", "length_increase=0.05 temperature=100", "length_increase=0.1 temperature=100",],
                         "strain",
                         "stress")


if __name__ == "__main__":
    # milestone4()
    milestone5_6()
    # milestone7()
    # milestone8()
    # milestone9()
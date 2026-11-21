
import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import legval


def legendre_polynomial(order, value):
    """Evaluate the Legendre polynomial P_order(value)."""
    coefficients = np.zeros(order + 1)
    coefficients[order] = 1.0
    return legval(value, coefficients)


def legendre_smoothing(y_data, olg=1, wnd=3):
    """
    Legendre polynomial smoothing.

    olg: number of Legendre polynomial terms
    wnd: half-width of the smoothing window
    """

    y_data = np.asarray(y_data, dtype=float)
    dt_sz = len(y_data)

    if wnd < 1:
        raise ValueError("wnd must be at least 1")

    if olg < 1:
        raise ValueError("olg must be at least 1")

    idm = 2 * wnd + 1

    if dt_sz < idm:
        raise ValueError(
            f"At least {idm} data points are required."
        )

    print("Maximum Legendre polynomial order (<olg):", olg - 1)
    print("Smoothing window (2*N+1):", idm)

    result = np.zeros(dt_sz)

    # Simpson integration weights
    weights = np.ones(idm)

    for j in range(1, idm - 1):
        weights[j] = 4 if j % 2 == 1 else 2

    # Coordinates within the smoothing window
    window_x = np.linspace(-1.0, 1.0, idm)

    # Calculate the Legendre expansion
    def calculate(window_y, dlt):

        smoothed_value = 0.0

        for order in range(olg):

            sm = 0.0

            for j in range(idm):
                sm += (
                    weights[j]
                    * legendre_polynomial(order, window_x[j])
                    * window_y[j]
                )

            sm *= (2 * order + 1) / (3.0 * (idm - 1))

            smoothed_value += (
                sm * legendre_polynomial(order, dlt)
            )

        return smoothed_value

    # -------------------------------------------------
    # Beginning of the data: first smoothing window
    # -------------------------------------------------

    first_window = y_data[:idm]

    for ic in range(wnd):

        dlt = (ic - wnd) / float(wnd)

        result[ic] = calculate(first_window, dlt)

    # -------------------------------------------------
    # Central part of the data
    # -------------------------------------------------

    for ic in range(wnd, dt_sz - wnd):

        window_y = y_data[ic - wnd:ic + wnd + 1]

        result[ic] = calculate(window_y, 0.0)

    # -------------------------------------------------
    # End of the data: last smoothing window
    # -------------------------------------------------

    last_window = y_data[dt_sz - idm:]

    for ic in range(dt_sz - wnd, dt_sz):

        dlt = (ic - dt_sz + 1 + wnd) / float(wnd)

        result[ic] = calculate(last_window, dlt)

    return result


def LgndrSld(filename="trial.dta", olg=5, wnd=7):

    # Read input data
    try:
        data = np.loadtxt(filename, ndmin=2)

    except OSError:
        print("Cannot open file:", filename)
        return

    except ValueError as error:
        print("Error reading data:", error)
        return

    if data.shape[1] < 2:
        print("Input file must contain at least two columns.")
        return

    # Original C++ format:
    # First column = yv2
    # Second column = xx

    y_data = data[:, 0]
    x_data = data[:, 1]

    if len(x_data) == 0:
        print("Input file contains no data.")
        return

    # Perform Legendre smoothing
    try:
        result = legendre_smoothing(y_data, olg, wnd)

    except ValueError as error:
        print("Smoothing error:", error)
        return

    # -------------------------------------------------
    # Plot original and smoothed data
    # -------------------------------------------------

    plt.figure(figsize=(10, 7))

    plt.plot(
        x_data,
        y_data,
        color="blue",
        linewidth=3,
        marker="o",
        markersize=4,
        label="Original data"
    )

    plt.plot(
        x_data,
        result,
        color="red",
        linewidth=2,
        marker="o",
        markersize=4,
        label="Legendre smoothed data"
    )

    plt.title("Legendre Smoothing")
    plt.xlabel("X")
    plt.ylabel("Y")

    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


# -----------------------------------------------------
# Main program
# -----------------------------------------------------

if __name__ == "__main__":

    filename = "trial.dta"
    olg = 5
    wnd = 11

    if len(sys.argv) > 1:
        filename = sys.argv[1]

    if len(sys.argv) > 2:
        olg = int(sys.argv[2])

    if len(sys.argv) > 3:
        wnd = int(sys.argv[3])

    LgndrSld(filename, olg, wnd)














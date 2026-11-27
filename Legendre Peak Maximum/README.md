# Cubic peak maximum fit

In the limit a peak maximum is quadratic about the maximum.  B
ut it is more practical, a better approach, to fit to a cubic and solve for the maximum:  Assume the peak is asymmetric.  

<img width="567" height="409" alt="PeakMaxFit" src="https://github.com/user-attachments/assets/87fedfe8-f9f7-4041-a57f-afd074e8995f" />. 

[Download ROOT Code file](https://github.com/user-attachments/files/35517078/LgndrPkMx.txt) for peak fit.  

Uses data file file from Legendre Smooth.   

import sys
import numpy as np
import matplotlib.pyplot as plt


def LgndrMax(filename="data.txt", iow=1, ipx=143):

    # -------------------------------------------------
    # Read input data
    # File format: yv2  xx
    # -------------------------------------------------

    try:
        data = np.loadtxt(filename)
    except OSError:
        print("Cannot open file:", filename)
        return

    if data.ndim != 2 or data.shape[1] < 2:
        print("Invalid input file format.")
        return

    yDt = data[:, 0].copy()
    x = data[:, 1].copy()

    dtSz = len(yDt)

    # -------------------------------------------------
    # Define peak window
    # -------------------------------------------------

    ipkSz = 5 + 2 * iow
    wnd = (ipkSz - 1) // 2

    start = ipx - 2 - iow
    end = start + ipkSz

    if start < 0 or end > dtSz:
        print("Peak window is outside the data range.")
        return

    xPk = x[start:end].copy()
    yPk = yDt[start:end].copy()

    dw = (xPk[-1] - xPk[0]) / float(ipkSz - 1)

    print(
        "step width", dw,
        "Peak window", ipkSz,
        "at", ipx,
        "file size", dtSz
    )

    # -------------------------------------------------
    # Calculate Legendre coefficients
    # Orders 0, 1, 2, 3
    # -------------------------------------------------

    lgr = np.zeros(4)

    for il in range(4):

        dlt = 0.0
        Sm = 0.0

        # First endpoint
        dlx = -1.0
        Sm += yPk[0] * np.polynomial.legendre.legval(
            dlx, [0] * il + [1]
        )

        # Last endpoint
        dlx = 1.0
        Sm += yPk[-1] * np.polynomial.legendre.legval(
            dlx, [0] * il + [1]
        )

        isW = abs(1)

        # Interior points
        for ica in range(1, ipkSz - 1):

            dlx = (-wnd + ica) / float(wnd)

            leg = np.polynomial.legendre.legval(
                dlx, [0] * il + [1]
            )

            Sm += 2 * (isW + 1) * leg * yPk[ica]

            isW -= 1
            isW = abs(isW)

        Sm *= (2 * il + 1) / (3.0 * float(ipkSz - 1))

        lgr[il] = Sm

        print(il, "lgndr cff", lgr[il])

    # -------------------------------------------------
    # Reconstruct fitted peak and calculate residuals
    # -------------------------------------------------

    Sum = 0.0
    yFit = np.zeros(ipkSz)

    for ica in range(ipkSz):

        dlx = (ica - wnd) / float(wnd)

        Sm = yPk[ica]

        # Legendre polynomial expansion
        yFit[ica] = (
            lgr[0] * np.polynomial.legendre.legval(dlx, [1])
            + lgr[1] * np.polynomial.legendre.legval(dlx, [0, 1])
            + lgr[2] * np.polynomial.legendre.legval(dlx, [0, 0, 1])
            + lgr[3] * np.polynomial.legendre.legval(dlx, [0, 0, 0, 1])
        )

        Sm -= yFit[ica]
        Sum += Sm * Sm

    Sum /= float(ipkSz - 5)

    # -------------------------------------------------
    # Calculate errors
    # -------------------------------------------------

    erry = np.sqrt(Sum)

    cb = lgr[3] * 2.5
    qd = lgr[2] * 1.5
    ln = lgr[1] - lgr[3] * 1.5

    if qd == 0:
        print("Error: qd is zero; cannot calculate errx.")
        return

    errx = np.sqrt(abs(erry / qd))

    # -------------------------------------------------
    # Determine maximum position
    # -------------------------------------------------

    discriminant = qd * qd - 3.0 * cb * ln

    if discriminant < 0:
        print("Error: no real peak maximum found.")
        return

    if cb == 0:
        print("Error: cubic coefficient is zero.")
        return

    sqrt_disc = np.sqrt(discriminant)

    xmx = (-qd - sqrt_disc) / (3.0 * cb)

    if abs(xmx) > 1.0:
        xmx = (-qd + sqrt_disc) / (3.0 * cb)

    # Calculate maximum amplitude
    ymx = (
        lgr[0] * np.polynomial.legendre.legval(xmx, [1])
        + lgr[1] * np.polynomial.legendre.legval(xmx, [0, 1])
        + lgr[2] * np.polynomial.legendre.legval(xmx, [0, 0, 1])
        + lgr[3] * np.polynomial.legendre.legval(xmx, [0, 0, 0, 1])
    )

    # Convert normalized coordinate to original x coordinate
    xmx = xPk[wnd] + xmx * wnd * dw

    # -------------------------------------------------
    # Print results
    # -------------------------------------------------

    print("xmx&xerr", xmx, errx)
    print("ymx&yerr", ymx, erry)

    # -------------------------------------------------
    # Plot original data and Legendre fit
    # -------------------------------------------------

    plt.figure(figsize=(8, 6))

    # Original data
    plt.plot(
        x, yDt,
        color="blue",
        linewidth=3,
        marker="o",
        markersize=4,
        label="Original data"
    )

    # Fitted peak
    plt.plot(
        xPk, yFit,
        color="red",
        linewidth=2,
        marker="o",
        markersize=4,
        label="Legendre cubic fit"
    )

    # Mark calculated maximum
    plt.plot(
        xmx, ymx,
        marker="*",
        color="green",
        markersize=12,
        label="Peak maximum"
    )

    plt.title("Legendre Peak Maximum Fit (Cubic)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()


# -----------------------------------------------------
# Main program
# -----------------------------------------------------

if __name__ == "__main__":

    filename = "data.txt"
    iow = 1
    ipx = 143

    if len(sys.argv) > 1:
        filename = sys.argv[1]

    if len(sys.argv) > 2:
        iow = int(sys.argv[2])

    if len(sys.argv) > 3:
        ipx = int(sys.argv[3])

    LgndrMax(filename, iow, ipx)



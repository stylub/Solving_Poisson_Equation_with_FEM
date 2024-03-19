package solution;

import com.github.sh0nk.matplotlib4j.Plot;
import com.github.sh0nk.matplotlib4j.PythonExecutionException;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class Main {
    private static final double G = 6.674e-11;

    public static int findIntervalId(double point, double h, int n) {
        for (int intervalID = 1; intervalID <= n; intervalID++) {
            double intervalBegin = h * (intervalID - 1), intervalEnd = h * intervalID;
            if (intervalBegin <= point && point <= intervalEnd) {
                return intervalID;
            }
        }
        return 1;
    }

    public static void main(String[] args) throws IOException, PythonExecutionException {
        int n = 200, begin = 0, end = 3, pointsOfInterest = 500;
        double h = (end - begin) / (double) n;

        // Matrix and Vector Initialization
        RealMatrix B = initializeMatrix(n, h);
        RealVector L = initializeVector(n, h);

        // Linear System Solution
        DecompositionSolver solver = new LUDecomposition(B).getSolver();
        RealVector solution = solver.solve(L);

        // Polynomial Function Initialization
        PolynomialFunction uHatted = new PolynomialFunction(new double[]{5, -1 / 3.0});

        List<Double> x = new ArrayList<>();
        List<Double> y = new ArrayList<>();

        // Data Generation Loop
        for (int counter = 0; counter < pointsOfInterest; counter++) {
            // Calculate point of interest
            double point = calculatePointOfInterest(counter, pointsOfInterest);

            // Identify the interval
            int intervalID = findIntervalId(point, h, n);

            // Evaluate Piecewise Function
            double phi = evaluatePiecewiseFunction(intervalID, point, solution, uHatted, h, n);

            x.add(point);
            y.add(phi);
        }

        Plot plt = Plot.create();
        plt.plot().add(x, y).label("MES");
        plt.legend().loc("upper right");
        plt.title(" Poisson equation for gravity in range [0 , 3] ");
        plt.show();
    }

    private static RealMatrix initializeMatrix(int n, double h) {
        double[][] matrixB = new double[n - 1][n - 1];

        Function function0 = new Function(h);

        for (int i = 0; i < n - 1; i++) {
            for (int j = 0; j < n - 1; j++) {
                matrixB[i][j] = computeMatrixElement(i, j, function0, h);
            }
        }

        return new Array2DRowRealMatrix(matrixB);
    }

    private static double computeMatrixElement(int i, int j, Function function0, double h) {
        if (i == j) {
            function0.setDerivativeFormula(i, j);
            function0.multiplyFormula(function0.getFunctionCoefficient());
            return 2 * (-1) * GaussianQuadrature.integrate(h * i, h * (i + 1), function0);
        } else if (i == j + 1 || j == i + 1) {
            Function derivativeProduct = new Function(new PolynomialFunction(new double[]{-1 / (h * h)}));
            return (-1) * GaussianQuadrature.integrate(h, h + h, derivativeProduct);
        } else {
            return 0;
        }
    }

    private static RealVector initializeVector(int n, double h) {
        double[] matrixL = new double[n - 1];

        Function function0 = new Function(h);
        Function function1 = new Function(h);

        for (int i = 0; i < n - 1; i++) {
            double beginCur = h * i, middleCur = h * (i + 1), endCur = h * (i + 2);

            matrixL[i] = computeLoadVectorElement(beginCur, middleCur, endCur, i, function0, function1, h);
        }

        return new ArrayRealVector(matrixL);
    }

    private static double computeLoadVectorElement(double beginCur, double middleCur, double endCur, int i,
                                                   Function function0, Function function1, double h) {
        double result = 0;

        if (1 <= middleCur && beginCur <= 2) {
            function0.setFormula(i + 1, i + 1);
            result += GaussianQuadrature.integrate(Math.max(1, beginCur), Math.min(2, middleCur), function0);
        }

        if (1 <= endCur && middleCur <= 2) {
            function0.setFormula(i + 1, i + 2);
            result += GaussianQuadrature.integrate(Math.max(1, middleCur), Math.min(2, endCur), function0);
        }

        function0.setDerivativeFormula(i, i);
        function0.multiplyFormula(-1 / 3.0);
        result -= (-1) * GaussianQuadrature.integrate(h * i, h * (i + 1), function0);

        function0.setDerivativeFormula(i, i + 1);
        function0.multiplyFormula(-1 / 3.0);
        result -= (-1) * GaussianQuadrature.integrate(h * (i + 1), h * (i + 2), function0);

        result *= Math.PI * 4 * G;

        return result;
    }

    private static double calculatePointOfInterest(int counter, int pointsOfInterest) {
        double increment = 3.0 / pointsOfInterest;
        return increment * counter;
    }

    private static double evaluatePiecewiseFunction(int intervalID, double point, RealVector solution,
                                                    PolynomialFunction uHatted, double h, int n) {
        Function function0 = new Function(h);
        Function function1 = new Function(h);

        double phi;
        if (intervalID == 1) {
            function0.setFormula(1, 1);
            phi = solution.getEntry(0) * function0.getFunction().value(point) + uHatted.value(point);
        } else if (intervalID == n) {
            function0.setFormula(n - 1, n);
            phi = solution.getEntry(n - 2) * function0.getFunction().value(point) + uHatted.value(point);
        } else {
            function0.setFormula(intervalID - 1, intervalID);
            function1.setFormula(intervalID, intervalID);
            phi = solution.getEntry(intervalID - 2) * function0.getFunction().value(point);
            phi += solution.getEntry(intervalID - 1) * function1.getFunction().value(point);
            phi += uHatted.value(point);
        }

        return phi;
    }
}



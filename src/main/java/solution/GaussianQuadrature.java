package solution;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;

import java.util.List;

public class GaussianQuadrature {
    private final int n = 5;
    private static final List<Double> weight = List.of(
            128.0 /225.0,
            (322 + 13 * Math.sqrt(70)) / 900,
            (322 + 13 * Math.sqrt(70)) / 900,
            (322 - 13 * Math.sqrt(70)) / 900,
            (322 - 13 * Math.sqrt(70)) / 900
    );
    private static final List<Double> points = List.of(
            0.0,
            1.0 / 3 * Math.sqrt(5 - 2 * Math.sqrt(10 / 7.0)),
            -1.0 / 3 * Math.sqrt(5 - 2 * Math.sqrt(10 / 7.0)),
            1.0 / 3 * Math.sqrt(5 + 2 * Math.sqrt(10 / 7.0)),
            -1.0 / 3 * Math.sqrt(5 + 2 * Math.sqrt(10 / 7.0))
    );
    static public double integrate(double begin, double end, Function function) {
        PolynomialFunction f = function.getFunction();
        double r = 0;
        int idx = 0;
        for (double point : points) {
            r += weight.get(idx) * f.value(((end - begin) / 2.0) * point + ((end + begin) / 2.0));
            idx++;
        }
        return r * ((end - begin) / 2.0);
    }
}
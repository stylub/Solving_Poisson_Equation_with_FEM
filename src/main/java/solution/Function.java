package solution;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;

interface FunctionFormula {
    double[] getCoefficients();
}

class ConstantFunction implements FunctionFormula {
    private final double coefficient;

    ConstantFunction(double coefficient) {
        this.coefficient = coefficient;
    }

    @Override
    public double[] getCoefficients() {
        return new double[]{coefficient};
    }
}

class LinearFunction implements FunctionFormula {
    private final double a;
    private final double b;

    LinearFunction(double a, double b) {
        this.a = a;
        this.b = b;
    }

    @Override
    public double[] getCoefficients() {
        return new double[]{a, b};
    }
}

public class Function {
    private final double h;
    private FunctionFormula functionFormula;

    public Function(double h) {
        this.h = h;
        initializeConstantFunction();
    }

    public Function(PolynomialFunction f) {
        this.h = 0;
        this.functionFormula = new ConstantFunction(f.getCoefficients()[0]);
    }

    public void setDerivativeFormula(int functionIntervalId, int wantedIntervalId) {
        FunctionFormula formula = getDerivativeFormula(functionIntervalId, wantedIntervalId);
        updateFunctionFormula(formula);
    }

    public void setFormula(int functionIntervalId, int wantedIntervalId) {
        FunctionFormula formula = getGeneralFormula(functionIntervalId, wantedIntervalId);
        updateFunctionFormula(formula);
    }

    public double getFunctionCoefficient() {
        return functionFormula.getCoefficients()[0];
    }

    public PolynomialFunction getFunction() {
        double[] coefficients = functionFormula.getCoefficients();
        return new PolynomialFunction(coefficients);
    }

    public void multiplyFormula(double number) {
        FunctionFormula multipliedFormula = new ConstantFunction(getFunctionCoefficient() * number);
        updateFunctionFormula(multipliedFormula);
    }

    private void initializeConstantFunction() {
        functionFormula = new ConstantFunction(0);
    }

    private void updateFunctionFormula(FunctionFormula formula) {
        functionFormula = formula;
    }

    private FunctionFormula getDerivativeFormula(int functionIntervalId, int wantedIntervalId) {
        if (functionIntervalId == wantedIntervalId) {
            return new ConstantFunction(1 / h);
        } else if (functionIntervalId + 1 == wantedIntervalId) {
            return new ConstantFunction(-1 / h);
        } else {
            return new ConstantFunction(0);
        }
    }

    private FunctionFormula getGeneralFormula(int functionIntervalId, int wantedIntervalId) {
        if (functionIntervalId == wantedIntervalId) {
            return new LinearFunction(1 - functionIntervalId, 1 / h);
        } else if (functionIntervalId + 1 == wantedIntervalId) {
            return new LinearFunction(1 + functionIntervalId, -1 / h);
        } else {
            return new ConstantFunction(0);
        }
    }
}

import java.io.FileReader;
import java.io.IOException;
import java.math.BigInteger;
import java.util.*;
import org.json.simple.*;
import org.json.simple.parser.*;

public class SecretSharingWithErrors {
    public static void main(String[] args) {
        // Read the JSON file
        JSONParser parser = new JSONParser();
        try {
            Object obj = parser.parse(new FileReader("input2.json"));
            JSONObject jsonObject = (JSONObject) obj;

            // Get n and k from "keys"
            JSONObject keys = (JSONObject) jsonObject.get("keys");
            int n = Integer.parseInt(keys.get("n").toString());
            int k = Integer.parseInt(keys.get("k").toString());

            // Collect x and y values
            List<BigInteger> xValues = new ArrayList<>();
            List<BigInteger> yValues = new ArrayList<>();

            // Iterate over the entries
            for (Object keyObj : jsonObject.keySet()) {
                String key = (String) keyObj;
                if (key.equals("keys")) {
                    continue;
                }
                // x is the key converted to BigInteger
                BigInteger x = new BigInteger(key);

                JSONObject entry = (JSONObject) jsonObject.get(key);
                int base = Integer.parseInt(entry.get("base").toString());
                String valueStr = entry.get("value").toString();

                // Decode y value
                BigInteger y = new BigInteger(valueStr, base);

                xValues.add(x);
                yValues.add(y);
            }

            // Ensure that n >= k + 2e
            int e = 1; // Maximum number of errors to correct
            if (xValues.size() < k + 2 * e) {
                System.err.println("Not enough data points to correct " + e + " errors.");
                return;
            }

            // Solve for P(x) and E(x)
            PolynomialResult result = berlekampWelch(xValues, yValues, k - 1, e);

            if (result == null) {
                System.err.println("Could not find a valid polynomial consistent with the data.");
                return;
            }

            // Extract the constant term c
            Rational c = result.getConstantTerm();

            if (c == null) {
                System.err.println("Unable to compute constant term c due to E(t) being zero at evaluation point.");
                return;
            }

            // Simplify c and print numerator if denominator is 1
            if (c.denominator.equals(BigInteger.ONE)) {
                System.out.println(c.numerator);
            } else {
                System.out.println("The constant term c is a fraction: " + c);
            }

            // Identify and print the wrong data points
            List<BigInteger> wrongDataPoints = result.findWrongDataPoints(xValues, yValues);
            if (!wrongDataPoints.isEmpty()) {
                System.out.println("Wrong data points detected at x values: " + wrongDataPoints);
            } else {
                System.out.println("No wrong data points detected.");
            }

        } catch (IOException e) {
            e.printStackTrace();
        } catch (ParseException e) {
            e.printStackTrace();
        }
    }

    private static PolynomialResult berlekampWelch(List<BigInteger> xValues, List<BigInteger> yValues, int k, int e) {
        int n = xValues.size();
    
        // Degrees of P(x) and E(x)
        int pDegree = k + e - 1; // Corrected degree of P(x)
        int eDegree = e;
    
        // Total number of unknowns (number of coefficients)
        int numUnknowns = (pDegree + 1) + (eDegree + 1); // Corrected calculation
    
        // Build the matrix for the linear system
        Rational[][] matrix = new Rational[n][numUnknowns];
    
        for (int i = 0; i < n; i++) {
            BigInteger xi = xValues.get(i);
            BigInteger yi = yValues.get(i);
    
            // Construct the row for the equation P(xi) = yi * E(xi)
            Rational[] row = new Rational[numUnknowns];
    
            // P(xi) coefficients (indices 0 to pDegree)
            for (int j = 0; j <= pDegree; j++) {
                row[j] = new Rational(xi.pow(j), BigInteger.ONE);
            }
    
            // -yi * E(xi) coefficients (indices pDegree + 1 to numUnknowns - 1)
            for (int j = 0; j <= eDegree; j++) {
                Rational term = new Rational(xi.pow(j).multiply(yi).negate(), BigInteger.ONE);
                row[pDegree + 1 + j] = term;
            }
    
            matrix[i] = row;
        }
    
        // Solve the linear system
        Rational[] solution = solveLinearSystem(matrix);
    
        if (solution == null) {
            return null;
        }
    
        // Extract P(x) and E(x) coefficients
        Rational[] pCoeffs = Arrays.copyOfRange(solution, 0, pDegree + 1);
        Rational[] eCoeffs = Arrays.copyOfRange(solution, pDegree + 1, solution.length);
    
        return new PolynomialResult(pCoeffs, eCoeffs);
    }
    
    private static Rational[] solveLinearSystem(Rational[][] matrix) {
        int n = matrix.length;
        int m = matrix[0].length - 1; // Exclude the augmented column

        // Perform Gaussian elimination with partial pivoting
        for (int col = 0; col < m; col++) {
            // Find the pivot row
            int pivotRow = -1;
            Rational maxVal = new Rational(BigInteger.ZERO, BigInteger.ONE);
            for (int row = col; row < n; row++) {
                Rational currentVal = matrix[row][col];
                if (currentVal.numerator.abs().compareTo(maxVal.numerator.abs()) > 0) {
                    maxVal = currentVal;
                    pivotRow = row;
                }
            }

            if (pivotRow == -1 || matrix[pivotRow][col].numerator.equals(BigInteger.ZERO)) {
                continue; // Cannot proceed with zero pivot
            }

            // Swap rows if necessary
            if (pivotRow != col) {
                Rational[] temp = matrix[col];
                matrix[col] = matrix[pivotRow];
                matrix[pivotRow] = temp;
            }

            // Normalize the pivot row
            Rational pivot = matrix[col][col];
            for (int j = col; j < matrix[col].length; j++) {
                matrix[col][j] = matrix[col][j].divide(pivot);
            }

            // Eliminate the current column entries in other rows
            for (int row = 0; row < n; row++) {
                if (row != col && !matrix[row][col].numerator.equals(BigInteger.ZERO)) {
                    Rational factor = matrix[row][col];
                    for (int j = col; j < matrix[row].length; j++) {
                        matrix[row][j] = matrix[row][j].subtract(matrix[col][j].multiply(factor));
                    }
                }
            }
        }

        // Extract the solution
        Rational[] solution = new Rational[m];
        for (int i = 0; i < m; i++) {
            solution[i] = matrix[i][matrix[i].length - 1];
        }

        return solution;
    }

    // Make the PolynomialResult class static
    static class PolynomialResult {
        public Rational[] pCoeffs;
        public Rational[] eCoeffs;

        public PolynomialResult(Rational[] pCoeffs, Rational[] eCoeffs) {
            this.pCoeffs = pCoeffs;
            this.eCoeffs = eCoeffs;
        }

        public Rational getConstantTerm() {
            // Choose x = 1 as the evaluation point
            BigInteger t = BigInteger.ONE;

            Rational pAtT = evaluatePolynomial(pCoeffs, t);
            Rational eAtT = evaluatePolynomial(eCoeffs, t);

            if (eAtT.numerator.equals(BigInteger.ZERO)) {
                System.err.println("E(" + t + ") is zero. Cannot compute constant term c at x = " + t);
                return null;
            }

            return pAtT.divide(eAtT);
        }

        public List<BigInteger> findWrongDataPoints(List<BigInteger> xValues, List<BigInteger> yValues) {
            List<BigInteger> wrongDataPoints = new ArrayList<>();

            for (int i = 0; i < xValues.size(); i++) {
                BigInteger xi = xValues.get(i);
                BigInteger yi = yValues.get(i);

                Rational pAtXi = evaluatePolynomial(pCoeffs, xi);
                Rational eAtXi = evaluatePolynomial(eCoeffs, xi);

                if (eAtXi.numerator.equals(BigInteger.ZERO)) {
                    // Avoid division by zero
                    wrongDataPoints.add(xi);
                    continue;
                }

                Rational expected = pAtXi.divide(eAtXi);

                if (!expected.equals(new Rational(yi, BigInteger.ONE))) {
                    wrongDataPoints.add(xi);
                }
            }

            return wrongDataPoints;
        }

        private Rational evaluatePolynomial(Rational[] coeffs, BigInteger x) {
            Rational result = new Rational(BigInteger.ZERO, BigInteger.ONE);
            Rational xPower = new Rational(BigInteger.ONE, BigInteger.ONE);

            for (Rational coeff : coeffs) {
                result = result.add(coeff.multiply(xPower));
                xPower = xPower.multiply(new Rational(x, BigInteger.ONE));
            }

            return result;
        }
    }

    // Make the Rational class static
    static class Rational {
        public BigInteger numerator;
        public BigInteger denominator;

        public Rational(BigInteger numerator, BigInteger denominator) {
            // Ensure denominator is positive
            if (denominator.signum() == -1) {
                numerator = numerator.negate();
                denominator = denominator.negate();
            }
            if (denominator.equals(BigInteger.ZERO)) {
                throw new ArithmeticException("Denominator cannot be zero.");
            }
            this.numerator = numerator;
            this.denominator = denominator;
            reduce();
        }

        public Rational add(Rational other) {
            BigInteger newNumerator = this.numerator.multiply(other.denominator).add(other.numerator.multiply(this.denominator));
            BigInteger newDenominator = this.denominator.multiply(other.denominator);
            return new Rational(newNumerator, newDenominator);
        }

        public Rational subtract(Rational other) {
            BigInteger newNumerator = this.numerator.multiply(other.denominator).subtract(other.numerator.multiply(this.denominator));
            BigInteger newDenominator = this.denominator.multiply(other.denominator);
            return new Rational(newNumerator, newDenominator);
        }

        public Rational multiply(Rational other) {
            BigInteger newNumerator = this.numerator.multiply(other.numerator);
            BigInteger newDenominator = this.denominator.multiply(other.denominator);
            return new Rational(newNumerator, newDenominator);
        }

        public Rational divide(Rational other) {
            if (other.numerator.equals(BigInteger.ZERO)) {
                throw new ArithmeticException("Cannot divide by zero.");
            }
            BigInteger newNumerator = this.numerator.multiply(other.denominator);
            BigInteger newDenominator = this.denominator.multiply(other.numerator);
            if (newDenominator.signum() == -1) {
                newNumerator = newNumerator.negate();
                newDenominator = newDenominator.negate();
            }
            return new Rational(newNumerator, newDenominator);
        }

        public void reduce() {
            BigInteger gcd = numerator.gcd(denominator);
            numerator = numerator.divide(gcd);
            denominator = denominator.divide(gcd);
        }

        @Override
        public boolean equals(Object obj) {
            if (obj instanceof Rational) {
                Rational other = (Rational) obj;
                return this.numerator.equals(other.numerator) && this.denominator.equals(other.denominator);
            }
            return false;
        }

        @Override
        public String toString() {
            return numerator + "/" + denominator;
        }
    }
}

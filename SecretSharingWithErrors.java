import java.io.FileReader;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.MathContext;
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
            int k = Integer.parseInt(keys.get("k").toString()); // k = degree + 1

            // Collect x and y values
            List<BigDecimal> xValues = new ArrayList<>();
            List<BigDecimal> yValues = new ArrayList<>();

            // Iterate over the entries
            for (Object keyObj : jsonObject.keySet()) {
                String key = (String) keyObj;
                if (key.equals("keys")) {
                    continue;
                }
                // x is the key converted to BigDecimal
                BigDecimal x = new BigDecimal(key);

                JSONObject entry = (JSONObject) jsonObject.get(key);
                int base = Integer.parseInt(entry.get("base").toString());
                String valueStr = entry.get("value").toString();

                // Decode y value
                BigInteger yInt = new BigInteger(valueStr, base);
                BigDecimal y = new BigDecimal(yInt);

                xValues.add(x);
                yValues.add(y);
            }

            // Ensure that n >= k + e
            int e = 3; // Maximum number of errors to correct
            if (xValues.size() < k + e) {
                System.err.println("Not enough data points to correct " + e + " errors.");
                return;
            }

            // Identify wrong data points and reconstruct the polynomial
            Result result = identifyAndReconstruct(xValues, yValues, k - 1, e);

            if (result == null) {
                System.err.println("Could not find a valid polynomial consistent with the data.");
                return;
            }

            // Extract the constant term c
            BigDecimal c = result.coefficients[0];

            // Simplify c if it's an integer
            if (c.stripTrailingZeros().scale() <= 0) {
                System.out.println(c.toBigIntegerExact());
            } else {
                System.out.println("The constant term c is: " + c.toPlainString());
            }

            // Print the wrong data points
            if (!result.wrongDataPoints.isEmpty()) {
                System.out.println("Wrong data points detected at indices: " + result.wrongDataPoints);
            } else {
                System.out.println("No wrong data points detected.");
            }

        } catch (IOException | ParseException e) {
            e.printStackTrace();
        }
    }

    private static Result identifyAndReconstruct(List<BigDecimal> xValues, List<BigDecimal> yValues, int degree, int maxErrors) {
        int n = xValues.size();
        List<Integer> indices = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            indices.add(i);
        }

        // Set precision for BigDecimal calculations
        MathContext mc = new MathContext(50);

        // Initial fit with all data points
        BigDecimal[] coefficients = fitPolynomial(xValues, yValues, degree, indices, mc);
        if (coefficients == null) {
            return null;
        }

        // Iteratively remove outliers
        for (int iteration = 0; iteration < maxErrors; iteration++) {
            // Calculate residuals
            Map<Integer, BigDecimal> residuals = new HashMap<>();
            for (int idx : indices) {
                BigDecimal xi = xValues.get(idx);
                BigDecimal yi = yValues.get(idx);
                BigDecimal yiFit = evaluatePolynomial(coefficients, xi, mc);
                BigDecimal residual = yi.subtract(yiFit, mc).abs(mc);
                residuals.put(idx, residual);
            }

            // Find the index with the maximum residual
            int worstIdx = Collections.max(residuals.entrySet(), Map.Entry.comparingByValue()).getKey();
            BigDecimal maxResidual = residuals.get(worstIdx);

            // Threshold to decide if it's an outlier (you may need to adjust this)
            if (maxResidual.compareTo(BigDecimal.ZERO) > 0) {
                // Remove the outlier and refit
                indices.remove((Integer) worstIdx);
                coefficients = fitPolynomial(xValues, yValues, degree, indices, mc);
                if (coefficients == null) {
                    return null;
                }
            } else {
                // No significant residuals, break
                break;
            }
        }

        // Wrong data points are those not in the final indices list
        List<Integer> wrongDataPoints = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            if (!indices.contains(i)) {
                wrongDataPoints.add(i);
            }
        }

        return new Result(coefficients, wrongDataPoints);
    }

    private static BigDecimal[] fitPolynomial(List<BigDecimal> xValues, List<BigDecimal> yValues, int degree, List<Integer> indices, MathContext mc) {
        int n = indices.size();
        int m = degree + 1;

        if (n < m) {
            return null; // Not enough data points
        }

        // Construct matrices
        BigDecimal[][] A = new BigDecimal[n][m];
        BigDecimal[] b = new BigDecimal[n];

        for (int i = 0; i < n; i++) {
            int idx = indices.get(i);
            BigDecimal xi = xValues.get(idx);
            BigDecimal yi = yValues.get(idx);

            for (int j = 0; j < m; j++) {
                A[i][j] = xi.pow(j, mc);
            }
            b[i] = yi;
        }

        // Solve the least squares problem
        return leastSquares(A, b, mc);
    }

    private static BigDecimal[] leastSquares(BigDecimal[][] A, BigDecimal[] b, MathContext mc) {
        // Use normal equations: (A^T * A) * x = A^T * b
        int m = A[0].length;

        // Compute A^T * A
        BigDecimal[][] ATA = new BigDecimal[m][m];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                ATA[i][j] = BigDecimal.ZERO;
                for (int k = 0; k < A.length; k++) {
                    ATA[i][j] = ATA[i][j].add(A[k][i].multiply(A[k][j], mc), mc);
                }
            }
        }

        // Compute A^T * b
        BigDecimal[] ATb = new BigDecimal[m];
        for (int i = 0; i < m; i++) {
            ATb[i] = BigDecimal.ZERO;
            for (int k = 0; k < A.length; k++) {
                ATb[i] = ATb[i].add(A[k][i].multiply(b[k], mc), mc);
            }
        }

        // Solve ATA * x = ATb
        return solveLinearSystem(ATA, ATb, mc);
    }

    private static BigDecimal[] solveLinearSystem(BigDecimal[][] A, BigDecimal[] b, MathContext mc) {
        int n = A.length;
        BigDecimal[][] augmentedMatrix = new BigDecimal[n][n + 1];

        // Build augmented matrix
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, augmentedMatrix[i], 0, n);
            augmentedMatrix[i][n] = b[i];
        }

        // Perform Gaussian elimination with partial pivoting
        for (int col = 0; col < n; col++) {
            // Find the pivot row
            int pivotRow = col;
            BigDecimal maxVal = augmentedMatrix[col][col].abs(mc);
            for (int row = col + 1; row < n; row++) {
                BigDecimal currentVal = augmentedMatrix[row][col].abs(mc);
                if (currentVal.compareTo(maxVal) > 0) {
                    maxVal = currentVal;
                    pivotRow = row;
                }
            }

            if (augmentedMatrix[pivotRow][col].compareTo(BigDecimal.ZERO) == 0) {
                return null; // Singular matrix
            }

            // Swap rows if necessary
            if (pivotRow != col) {
                BigDecimal[] temp = augmentedMatrix[col];
                augmentedMatrix[col] = augmentedMatrix[pivotRow];
                augmentedMatrix[pivotRow] = temp;
            }

            // Normalize the pivot row
            BigDecimal pivot = augmentedMatrix[col][col];
            for (int j = col; j <= n; j++) {
                augmentedMatrix[col][j] = augmentedMatrix[col][j].divide(pivot, mc);
            }

            // Eliminate the current column entries in other rows
            for (int row = 0; row < n; row++) {
                if (row != col) {
                    BigDecimal factor = augmentedMatrix[row][col];
                    for (int j = col; j <= n; j++) {
                        augmentedMatrix[row][j] = augmentedMatrix[row][j].subtract(augmentedMatrix[col][j].multiply(factor, mc), mc);
                    }
                }
            }
        }

        // Extract the solution
        BigDecimal[] solution = new BigDecimal[n];
        for (int i = 0; i < n; i++) {
            solution[i] = augmentedMatrix[i][n];
        }

        return solution;
    }

    private static BigDecimal evaluatePolynomial(BigDecimal[] coeffs, BigDecimal x, MathContext mc) {
        BigDecimal result = BigDecimal.ZERO;
        BigDecimal xPower = BigDecimal.ONE;

        for (BigDecimal coeff : coeffs) {
            result = result.add(coeff.multiply(xPower, mc), mc);
            xPower = xPower.multiply(x, mc);
        }

        return result;
    }

    static class Result {
        public BigDecimal[] coefficients;
        public List<Integer> wrongDataPoints;

        public Result(BigDecimal[] coefficients, List<Integer> wrongDataPoints) {
            this.coefficients = coefficients;
            this.wrongDataPoints = wrongDataPoints;
        }
    }
}

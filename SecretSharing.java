import java.io.FileReader;
import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import org.json.simple.*;
import org.json.simple.parser.*;

public class SecretSharing {
    public static void main(String[] args) {
        // Read the JSON file
        JSONParser parser = new JSONParser();
        try {
            Object obj = parser.parse(new FileReader("input.json"));
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

            if (xValues.size() < k) {
                System.err.println("Not enough data points to reconstruct the polynomial.");
                return;
            }

            // Use the first k points
            xValues = xValues.subList(0, k);
            yValues = yValues.subList(0, k);

            // Calculate c using Lagrange interpolation at x = 0
            Rational c = lagrangeInterpolationAtZero(xValues, yValues);

            // Simplify c and print numerator if denominator is 1
            if (c.denominator.equals(BigInteger.ONE)) {
                System.out.println(c.numerator);
            } else {
                System.out.println("The constant term c is a fraction: " + c);
            }

        } catch (IOException e) {
            e.printStackTrace();
        } catch (ParseException e) {
            e.printStackTrace();
        }
    }

    private static Rational lagrangeInterpolationAtZero(List<BigInteger> xValues, List<BigInteger> yValues) {
        Rational c = new Rational(BigInteger.ZERO, BigInteger.ONE);
        int k = xValues.size();

        for (int i = 0; i < k; i++) {
            Rational term = new Rational(yValues.get(i), BigInteger.ONE);
            for (int j = 0; j < k; j++) {
                if (i != j) {
                    Rational fraction = new Rational(xValues.get(j).negate(), xValues.get(i).subtract(xValues.get(j)));
                    term = term.multiply(fraction);
                }
            }
            c = c.add(term);
        }

        c = c.reduce();
        return c;
    }
}

class Rational {
    public BigInteger numerator;
    public BigInteger denominator;

    public Rational(BigInteger numerator, BigInteger denominator) {
        // Ensure denominator is positive
        if (denominator.signum() == -1) {
            numerator = numerator.negate();
            denominator = denominator.negate();
        }
        this.numerator = numerator;
        this.denominator = denominator;
    }

    public Rational add(Rational other) {
        BigInteger newNumerator = this.numerator.multiply(other.denominator).add(other.numerator.multiply(this.denominator));
        BigInteger newDenominator = this.denominator.multiply(other.denominator);
        return new Rational(newNumerator, newDenominator).reduce();
    }

    public Rational multiply(Rational other) {
        BigInteger newNumerator = this.numerator.multiply(other.numerator);
        BigInteger newDenominator = this.denominator.multiply(other.denominator);
        return new Rational(newNumerator, newDenominator).reduce();
    }

    public Rational reduce() {
        BigInteger gcd = numerator.gcd(denominator);
        return new Rational(numerator.divide(gcd), denominator.divide(gcd));
    }

    @Override
    public String toString() {
        return numerator + "/" + denominator;
    }
}

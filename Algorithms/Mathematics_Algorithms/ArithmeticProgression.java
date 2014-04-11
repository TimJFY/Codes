package Mathematics_Algorithms;

// compute the sum of the first n positive natural numbers without using any of the following:
// if-else, while, for, foreach, switch-case, '? : ', '<'('<='), '>'('>=')
public class ArithmeticProgression {
	static A objA;
	static A objB;

	public static void main(String[] args) {
		ArithmeticProgression arithmeticProgression = new ArithmeticProgression();
		objA = arithmeticProgression.new A();
		objB = arithmeticProgression.new B();
		objB.object[0] = objA;
		objB.object[1] = objB;
		System.out.println("Sum = " + objB.sum(50));
	}

	class A {
		A object[] = new A[2];

		public int sum(int n) {
			return 0;
		}
	}

	class B extends A {
		public int sum(int n) {
			// In C++ : object[!!n]
			boolean b = (n != 0);
			return object[-("false".indexOf("" + b))].sum(n - 1) + n;
		}
	}
}

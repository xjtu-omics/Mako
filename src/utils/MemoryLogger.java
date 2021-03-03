package utils;


/**
 * This class is used to record the maximum memory usage of the algorithm during execution.
 */
public class MemoryLogger {
	
	// the only instance  of this class (this is the "singleton" design pattern)
	private static MemoryLogger instance = new MemoryLogger();

	// variable to store the maximum memory usage
	private double maxMemory = 0;
	
	/**
	 * Method to obtain the only instance of this class
	 * @return instance of MemoryLogger
	 */
	public static MemoryLogger getInstance(){
		return instance;
	}
	
	/**
	 * To get the maximum amount of memory used until now
	 * @return a double value indicating memory as megabytes
	 */
	public double getMaxMemory() {
		return maxMemory;
	}

	/**
	 * Reset the maximum amount of memory recorded.
	 */
	public void reset(){
		maxMemory = 0;
	}
	
	/**
	 * Check the current memory usage and record it if it is higher
	 * than the amount of memory previously recorded.
	 */
	public void checkMemory() {
		double currentMemory = (Runtime.getRuntime().totalMemory() -  Runtime.getRuntime().freeMemory())
				/ 1024d / 1024d;
		if (currentMemory > maxMemory) {
			maxMemory = currentMemory;
		}
	}
}

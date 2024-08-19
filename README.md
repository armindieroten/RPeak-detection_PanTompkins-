# ECG Real-Time Plotting

 **ECG_Real-Time_Plotting.ipynb** enables real-time ECG signal visualization by using parallel threads to simultaneously write to and read from a serial port. The signal, along with its detected peaks, is rendered in real time using the `pyqtgraph` library.

## How to Use:

1. Install the required dependencies by running:
   ```bash
   !pip install -r requirements.txt
2. In `ECG_Real-Time_Plotting.ipynb` modify the `serial.Serial()` parameters to match your configuration (e.g., COM port, baud rate, parity, etc.).

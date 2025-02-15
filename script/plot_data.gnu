# Set output format
set terminal pngcairo enhanced size 1200,800
set output 'temperature_kinetic_energy.png'

# Labels and title
set title "Temperature and Kinetic Energy Over Time"
set xlabel "Time (fs)"
set ylabel "Kinetic Energy"
set y2label "Temperature (K)"

# Configure axes
set y2tics     # Enable right y-axis ticks
set grid       # Enable grid

# Load the data
set datafile separator ","

# Plot the data with two axes
plot "../out/particules_data.csv" using 1:3 with lines title "Kinetic Energy" lw 2 lc rgb "blue", \
     "../out/particules_data.csv" using 1:2 axes x1y2 with lines title "Temperature (K)" lw 2 lc rgb "red"

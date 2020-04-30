The files are named as:

"open_ch1", "open_ch4", "open_ch16", "apt_ch1", "apt_ch4", "apt_ch16", 

the first part "open" or "apt" indicates the experiment environment. 
This is followed by the channel number: "ch1", "ch4" or "ch16". 

Each environment contains three experiments (with different number of channels) and in each experiment, four different trials are performed, each with a unique trajectory.
In the first trial, trial = 0, the person continuously walks around and the ground truth trajectory is unknown.
In trials 1-3, the person walks from one reference position to another and in each position, the person remains stationary for a few seconds. The ground truth trajectory is known for these trials.

The one channel system uses channel 26.
The four channel system uses channels 26,16,21,11
The sixteen channel system uses channels 26,22,18,14,25,21,17,13,24,20,16,12,23,19,15,11

The coordinates of the nodes are:

OPEN ENVIRONMENT NODE LOCATIONS:

    x-coor    y-coor
  
    2.0200         0
    3.4200         0
    4.8200         0
    6.2200         0
    7.6200         0
    9.5800    1.7900
    9.5800    2.8400
    9.5800    3.8900
    9.5800    4.9400
    9.5800    5.9900
    7.6200    7.8200
    6.2200    7.8200
    4.8200    7.8200
    3.4200    7.8200
    2.0200    7.8200
         0    5.9900
         0    4.9400
         0    3.8900
         0    2.8400
         0    1.7900


APARTMENT ENVIRONMENT NODE LOCATIONS:

    x-coor    y-coor

    1.9600    0.2700
    0.0800    2.8900
    0.0800    4.5400
    0.0800    6.5700
    3.1000    5.8200
    2.2600    3.9500
    3.1200    2.6100
    3.3800    5.9900
    5.5300    8.0000
    6.9000    5.3100
    6.5200    3.7000
    6.4000    0.0800
    5.0200    0.0800
    7.1800    6.2000
    8.3400    7.4900
   10.2400    6.9200
   10.2400    4.3800
    7.9000    3.9300
   10.2400    2.8000
    8.3400    0.2700

# The complexity of geometric scaling
This is the Python implementation of two geometric scaling algorithms: **MRA-based geometric scaling** and **Geometric Scaling via a feasibility test**.

## Setups
All code was developed and tested on a MacBook Pro with a 2GHz i5 processor and 16GB of RAM. The environment is as bellow:

- macOS Monterey
- Python 3.9.7 (Anaconda 4.12.0 64 bit)


## How to Run
Here are two examples:
```bash
python3 main.py --algorithm 1 --N 20 --visualization 1 --earlystopping 0 
```
```bash
python3 main.py --algorithm 2 --N 20 --visualization 1 --earlystopping 1 --molecular 4 --denominator 3
```

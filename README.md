# ClimateHPC
Final project for HPC4DS course held at UniTn
This is the structure of the files contained in the repo.
Each folder in final settings contains a bash file and a c file for the two kinds of files.
For the plotting directory you can find all the pictures for the 12 average velocity and sea surface elevation for the 12 month in high and low resolutions.

This files were designed to work for High performance computers of UNITN and trying to run these codes on any different computer would require reformalation of the bash files used to submit the jobs. 
```bash
ClimateHPC/
├── Final_settings
│   ├── elevation_data
│   │   ├── reading_ssh_early.c
│   │   ├── reading_ssh_early.sh
│   │   ├── reading_ssh_f.c
│   │   └── reading_ssh_f.sh
│   └── velocity_data
│       ├── reading_u_hybrid_v1.c
│       ├── reading_u_hybrid_v2.c
│       ├── reading_u_simple.c
│       ├── reading_u_spill_comm_v2.c
│       ├── run_serial.sh
│       ├── run.sh
│       └── serial.c
├── plotting
│   ├── final_plots
│   ├── low_res_plots
│   ├── map_plotter.py
│   └── tex_plotter.py
├── README.md
└── utils
    ├── tdata.csv
    └── utils.py
```

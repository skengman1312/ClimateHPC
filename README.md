# ClimateHPC
Final project for HPC4DS course held at UniTn
This is the structure of the files contained in the repo.
Each folder in final settings contains a bash file and a c file for the two kinds of files.
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

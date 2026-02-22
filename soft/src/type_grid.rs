use ndarray::{Array2};
use std::fs::File;
use std::io::{self, BufRead, BufReader};

use crate::constants::*;

pub struct Vec2 {
    pub(crate) x : f64,
    pub(crate) y : f64,
}

impl Vec2{
    fn new(x: f64, y: f64) -> Self {
    Self{x,y}
    }
}

pub(crate) enum GridType {
    CartMesh,
    GeneralQuads,
}

pub(crate) struct Grid {
    //Grid stuff and things
    //
    // Identified for what kind of grid this is
    // CURRENTLY ONLY CARTESEAN IS SUPPORTED
    pub(crate)itype : GridType
    // Contains all X, Y points of the grid
    // In i-major order. So:
    // for j=0->nj-1{ for ... i { x[i,j] = coords[n].x ; n++ } }
    pub(crate) coords : Vec<Vec2>,
    //    The convention is that arrays/dimensions are for numerical\conceptual
    //        ideas whereas the vector subtype is for physical vectors.
    //        I.E. that's the reason for this structure instead of a 2d array
    //        ... also this makes sure we have mem locality for xy points

    // The number of poins in each grid direction
    pub(crate) ni : usize,
    pub(crate) nj : usize,
    pub(crate) npoin : uzise,
    pub(crate) nelem : uzise,

    // Volume of cells (doing a cart grid for now 
    //               so just need one val suckers)
    pub(crate) vol : f64,

    // Face lengths 
    pub(crate) fleni : f64,
    pub(crate) flenj : f64,
    /// NOTE: for non-cart meshes we'll need face normal information
}

pub(crate) impl Grid {
    pub(crate) fn new(grid_file: File) -> Self {
        //
        // Load grid
        read_grid_file(grid_file)
    }

    fn read_grid_file(grid_file: File) -> Grid {
        // Read in grid size, type
        // Then coordinates (quads), or dims (cart)
        //
        // IMPORTANT: EVEN THOUGH WE READ IN GRID TYPE, THIS ROUTINE ONLY WORKS
        //              FOR CART GRIDS AND DOES NOT ADAPT TO OTHER GRID TYPES
        
        // Corners of grid
        let mut p0 = Vec2::new(0.0,0.0);
        let mut p1 = Vec2::new(0.0,0.0);

        // Open er up
        let reader = BufReader::new(grid_file);

        // Loop through all lines of file
        let lnum:usize = 0;
        for line in reader.lines() {
            lnum += 1 ;
            // Interpret line based on what number it is
            match lnum {
                2  => grid_type_index = line.as_usize(),
                4  => ni = line.as_f64(),
                5  => nj = line.as_f64(),
                7  => p0.x = line.as_f64(),
                8  => p0.y = line.as_f64(),
                10 => p1.x = line.as_f64(),
                11 => p1.y = line.as_f64(),
                _ => {}
            }
        }

        assert!(p1.x>p0.x and p1.y>p0.y, "donkey")
        //Now have all information we need to assemble a cartesean grid
        let npoin = ni * nj;
        let nelem = (ni-1) * (nj-1);

        // Face length
        let fleni = (p1.x-p0.x) / (ni-1);
        let flenj = (p1.y-p0.y) / (nj-1);

        // Cell Volume
        let vol = fleni * flenj;
        
        // Compute point coordinates
        let ptnum: usize = 0;
        let coords : Vec<Vec2> = [];
        for j in 0..nj {
            for i in 0..ni {
                let xpoint:f64 = i*p0.x + (1-i/(ni-1))*p1.x;
                let ypoint:f64 = j*p0.y + (1-j/(nj-1))*p1.y;
                coords.push(Vec2::new(xpoint,ypoint));
            }
        }

        Grid {ni,nj,npoin,nelem,fleni,flenj,vol,coords,
            itype : GridType::CartMesh}

    }
}

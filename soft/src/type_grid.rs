use ndarray::{Array2};
use std::fs::File;

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
    Cartesean_Mesh,
    General_Quads,
}

pub(crate) struct Grid {
    //Grid stuff and things
    //
    // Identified for what kind of grid this is
    // CURRENTLY ONLY CARTESEAN IS SUPPORTED
    pub(crate)itype : GridType
    // Contains all X, Y points of the grid
    // In j-major order. So:
    // for j=0->nj-1{ for ... i { x[i,j] = coords[n].x ; n++ } }
    pub(crate) coords : Array::<Vec2>,

    // The number of poins in each grid direction
    pub(crate) ni : usize,
    pub(crate) nj : usize,

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
        read_grid_file(grid_file);
        //
        //
        // Calculate volume and face lengths

    }

    fn read_grid_file(grid_file: File) {
        // Read in grid size, type
        // Then coordinates (quads), or dims (cart)

    }
}

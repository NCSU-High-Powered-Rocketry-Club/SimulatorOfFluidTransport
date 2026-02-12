use ndarray::{Array2};

use crate::constants::*;

pub(crate) struct Solution_Array {
    
    // Since array for holding the data
    data : Array2::<f64>,
}

impl SolutionArray {

    pub fn new(num_points : usize, num_variables : usize) -> Self {
        Self{data : Array2::<f64>::zeros((num_points, num_variables))} 
    }
    
    pub fn density(&self) -> ArrayView1<f64> {
        self.data.column(0)
    }
    pub fn density_mut(&mut self) -> ArrayViewMut1<f64> {
        self.data.column_mut(0)
    }
    pub fn density_owned(&self) -> Array1<f64> {
        self.data.column(0).to_owned()
    }

    pub fn momentum(&self) -> Array1::<f64> {
        self.column(1)
    }   
    
    pub fn energy(&self) -> Array1::<f64> {
        self.column(2)
    }
}

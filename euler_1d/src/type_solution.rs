use ndarray::{Array1,Array2,ArrayView1,ArrayViewMut1};

use crate::constants::*;

pub(crate) struct SolutionArray {
    
    // Since array for holding the data
    pub(crate) data : Array2::<f64>,
}

impl SolutionArray {

    pub fn new(num_points : usize, num_variables : usize) -> Self {
        Self{data : Array2::<f64>::zeros((num_points, num_variables))} 
    }
    
    pub fn column(&self, c : usize) -> ArrayView1<'_,f64> {
        self.data.column(c)
    }
    pub fn column_mut(&mut self, c : usize) -> ArrayViewMut1<'_,f64> {
        self.data.column_mut(c)
    }
    pub fn column_owned(&self, c : usize) -> Array1<f64> {
        self.data.column(c).to_owned()
    }

    pub fn row(&self, c : usize) -> ArrayView1<'_,f64> {
        self.data.row(c)
    }
    pub fn row_mut(&mut self, c : usize) -> ArrayViewMut1<'_,f64> {
        self.data.row_mut(c)
    }
    pub fn row_owned(&self, c : usize) -> Array1<f64> {
        self.data.row(c).to_owned()
    }
}

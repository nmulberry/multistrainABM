extern crate csv;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate rand;
extern crate rand_distr; 
use rand_distr::{Bernoulli, Distribution};
extern crate multiabm; // load library with model defn
use multiabm::*;
extern crate ode_event_solvers;
use ode_event_solvers::euler::*;
use ode_event_solvers::*;
use std::env;
//========================================================//
/*                 SETUP                                  */
//========================================================//
type State = DVector<f64>;

#[derive(Deserialize, Debug)]
struct SeroParameters {
  alpha: Vec<f64>,
}


impl SeroParameters {
  fn new() -> Self {
    SeroParameters {
      alpha: Vec::new(),
    }
  }

  fn push(&mut self, row: &csv::StringRecord) {
      self.alpha.push(row[0].parse().unwrap());
    }

    fn read_csv(filepath: &str) -> Self {// no header
      let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false).from_path(filepath).unwrap();
      let mut df = SeroParameters::new();
      for result in rdr.records() {
        df.push(&result.unwrap());
      }
      df
    }
}

#[derive(Deserialize, Debug)]
struct SCParameters {
  kappa_m: Vec<f64>,
}


impl SCParameters {
  fn new() -> Self {
    SCParameters {
      kappa_m: Vec::new(),
    }
  }

  fn push(&mut self, row: &csv::StringRecord) {
      self.kappa_m.push(row[0].parse().unwrap());
    }

    fn read_csv(filepath: &str) -> Self {// no header
      let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false).from_path(filepath).unwrap();
      let mut df = SCParameters::new();
      for result in rdr.records() {
        df.push(&result.unwrap());
      }
      df
    }
}

#[derive(Deserialize, Debug)]
struct AGParameters {
  kappa_o: Vec<f64>,
}


impl AGParameters {
  fn new() -> Self {
    AGParameters {
      kappa_o: Vec::new(),
    }
  }

  fn push(&mut self, row: &csv::StringRecord) {
      self.kappa_o.push(row[0].parse().unwrap());
    }

    fn read_csv(filepath: &str) -> Self {// no header
      let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false).from_path(filepath).unwrap();
      let mut df = AGParameters::new();
      for result in rdr.records() {
        df.push(&result.unwrap());
      }
      df
    }
}


#[derive(Deserialize, Debug)]
struct ConstParameters {
  t_max: f64,
  t_vax: f64,
  nhost: usize,
  nstrain: usize,
  beta: f64,
  theta: f64,
  eps_x: f64,
  eps_a: f64,
  rho: f64,
  t_g: f64,
  t_l: f64
}

impl ConstParameters {
  fn new() -> Self {
    ConstParameters {
      t_max: 0.0,
      t_vax: 0.0,
      nhost: 0,
      nstrain: 0,
      beta: 0.0,
      theta: 0.0,
      eps_x: 0.0,
      eps_a: 0.0,
      rho: 0.0,
      t_g: 0.0,
      t_l: 0.0, 
    }
  }

  fn read_csv(filepath: &str) -> Self {
    let mut rdr = csv::Reader::from_path(filepath).unwrap();
    let mut df = ConstParameters::new();
    let header = csv::StringRecord::from(vec![
        "t_max", "t_vax", "nhost", "nstrain", "beta", 
        "theta", "eps_x", "eps_a", "rho", "t_g", "t_l",
    ]);
    for result in rdr.records() {
      df = result.unwrap().deserialize(Some(&header)).unwrap()
    }
    df
  }
}

#[derive(Deserialize, Debug)]
struct StrainParameters {
  id_m: Vec<u64>,
  id_s: Vec<u64>,
  id_v: Vec<u64>,
  id_o: Vec<u64>,
  id_init: Vec<u64>,
}

impl StrainParameters {
  fn new() -> Self {
    StrainParameters {
      id_m: Vec::new(),
      id_s: Vec::new(),
      id_v: Vec::new(),
      id_o: Vec::new(),
      id_init: Vec::new(),
     }
  }

  fn push(&mut self, row: &csv::StringRecord) {
    self.id_m.push(row[0].parse().unwrap());
    self.id_s.push(row[1].parse().unwrap());
    self.id_o.push(row[2].parse().unwrap());
    self.id_v.push(row[3].parse().unwrap());
    self.id_init.push(row[4].parse().unwrap());
  }

  fn read_csv(filepath: &str) -> Self {
    let mut rdr = csv::Reader::from_path(filepath).unwrap();
    let mut df = StrainParameters::new();
    for result in rdr.records() {
      df.push(&result.unwrap());
    }
    df
  }
}

//========================================================//
/*                 MAIN                                  */
//========================================================//
fn main() {
    let strain_pars = StrainParameters::read_csv(&String::from("./strain_pars.csv"));
    let pars  = ConstParameters::read_csv(&String::from("./const_pars.csv"));
    let args: Vec<String> = env::args().collect();
    let sc_pars  = SCParameters::read_csv(&args[1]); //FNAME
    let sero_pars  = SeroParameters::read_csv(&args[2]); //FNAME
    let ag_pars  = AGParameters::read_csv(&args[3]); //FNAME

    // Initial state.   
    let mut y0 = State::from_vec(vec![0.0; pars.nhost*pars.nstrain*2]);
    // Choose hosts at random to be innoculated for each strain
    let mut i0 = 50_f64/(pars.nhost as f64); // high
    if i0 > 1.0 {i0 = 1.0;} // sample prob per host per strain
    let init_infect_dist = Bernoulli::new(i0).unwrap();

    for k in 0..pars.nhost {
      // pick hosts to be initially infected
      for j in 0..pars.nstrain {
        if strain_pars.id_init[j] == 1{
          let v = init_infect_dist.sample(&mut rand::thread_rng());
          if v {
            y0[k*pars.nstrain+j] = pars.rho;
          }
        }
      }
    }
    // compute kappa and alpha vecs (for each strain)
    let mut kappa_vec = vec![0.0; pars.nstrain];
    for j in 0..pars.nstrain {
      let s = strain_pars.id_m[j]-1;
      kappa_vec[j] += sc_pars.kappa_m[s as usize];
      let s = strain_pars.id_o[j]-1;
      kappa_vec[j] += ag_pars.kappa_o[s as usize];
    }
    let mut alpha_vec = vec![0.0; pars.nstrain];
    for j in 0..pars.nstrain {
      let s = strain_pars.id_s[j]-1;
      alpha_vec[j] = sero_pars.alpha[s as usize];
    }

    // Create structure to pass to DE solver
    let mut system = WhCompetition {
      nhost: pars.nhost, 
      nstrain: pars.nstrain,
      theta: pars.theta,
      eps_x: pars.eps_x,
      eps_a: pars.eps_a,
      rho: pars.rho,
      alpha: alpha_vec.clone(),
      kappa: kappa_vec.clone(),
      id_m: strain_pars.id_m,
      id_s: strain_pars.id_s,
      id_o: strain_pars.id_o,
      beta: pars.beta,
      t_g: pars.t_g,
      t_l: pars.t_l,
    }; 
   /* SOLVE SYSTEM */
    // pre-vaccination:
    let mut stepper = Euler::new(system.clone(), 0., y0, pars.t_vax, Vec::from([0.1, 1.0, 10.0]));
    let _res = stepper.integrate().unwrap(); //returns initial and final states
    if pars.t_vax < pars.t_max {
      // post-vaccination
      let x0 = stepper.x_out()[1].clone()+1.0;
      y0 = stepper.y_out()[1].clone();
      // add vaccine (growth rate to zero--> instable to bump up a)
      for j in 0..pars.nstrain {
        if strain_pars.id_v[j] == 1 {
          kappa_vec[j] = 0.5; // dies off eventually but not immediately
        } 
      }
      system.kappa = kappa_vec;
      stepper = Euler::new(system.clone(), x0, y0, pars.t_max, Vec::from([0.1, 1.0, 10.0]));
      let _res = stepper.integrate().unwrap();
    }
}
 

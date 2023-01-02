extern crate ode_event_solver;
use ode_event_solver::*;
extern crate rand;
extern crate rand_distr; 
use rand_distr::{Bernoulli, Poisson, Distribution};
use rand::distributions::WeightedIndex;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rand::seq::index::sample;

//========================================================//
/*                 MODEL DEFINITION                      */
//========================================================//
#[derive(Clone, Debug)]
pub struct WhCompetition {
  pub nhost: usize,
  pub nstrain: usize,
  pub tau: f64,
  pub theta: f64,
  pub eps_x: f64,
  pub eps_a: f64,
  pub rho: f64,
  pub alpha: Vec<f64>,
  pub kappa: Vec<f64>,
  pub id_m: Vec<u64>,
  pub id_s: Vec<u64>,
  pub id_r: Vec<u64>,
  pub p_treat: Bernoulli, 
  pub beta: f64,
  pub mu: f64,
  pub on_treat: Vec<bool>, 
}

type State = DVector<f64>;
type Time = f64;

// Within Host Model for arbitrary number of strains and hosts
impl ode_event_solver::System<State> for WhCompetition {
    fn ode(&self, _: Time, y: &State, dy: &mut State) {
      for k in 0..(self.nhost) {
        // For each strain in host:
        let mut m = vec![0.0; self.nstrain];
        let mut s = vec![0.0; self.nstrain];

        for i in 0..(self.nstrain) {
           let m_i = self.id_m[i];
           let s_i = self.id_s[i];
           let mut mm = 0.0;
           let mut ss = 0.0;
           for j in 0..(self.nstrain){
             if m_i==self.id_m[j] {
               mm += y[k*self.nstrain+j];
             }
             if s_i==self.id_s[j] {
               ss += y[k*self.nstrain+j];
             }
            }
           m[i] = mm;
           s[i] = ss;
       }

       for j in 0..(self.nstrain) {
          let ii = k*self.nstrain+j as usize;
          let jj = k*self.nstrain + self.nhost*self.nstrain+j as usize;
          // prevalence
          dy[ii] = y[ii]*(1.0-m[j])*self.kappa[j]- 
                    y[ii]*y[jj]*self.alpha[j]-
                    y[ii]*self.eps_x; 
          if self.on_treat[k] {
            if self.id_r[j] == 0 {
              dy[ii] = dy[ii] - y[ii]*self.tau;
            }        
          }

          // immunity
          dy[jj] = s[j]*s[j]/(self.rho*self.rho+s[j]*s[j])-
                     self.eps_a*y[jj]+
                     y[jj]*y[jj]/(self.theta*self.theta+
                                  y[jj]*y[jj]);
      }
    }
  }


  fn event(&mut self, _: Time, y: &State, dy: &mut State) {//abs freq
    // ynew = y + dy
    // Get carriage by strain. If density < rho, set to 0.
    let migrate_dist = Bernoulli::new(self.mu).unwrap(); // 1% per day
    let mut n_infected:Vec<f64> = vec![0.0; self.nstrain];
    for k in 0..self.nhost {
      // get carriage
      for j in 0..self.nstrain {
        let ii = k*self.nstrain+j; //density of strain j
        if y[ii] < self.rho {
          dy[ii] = -y[ii];
        } else {
          n_infected[j] += y[ii];  //total density of the strain
        }
      }
    }
    for j in 0..self.nstrain {
      if n_infected[j] > 0.0 {
        let transmit_dist = Poisson::new(self.beta*n_infected[j]).unwrap();
        let n_transmit = transmit_dist.sample(&mut rand::thread_rng());
        let recip = sample(&mut rand::thread_rng(), self.nhost, n_transmit as usize);
        for r in recip {
           self.on_treat[r] = self.p_treat.sample(&mut rand::thread_rng());
           dy[r*self.nstrain+j] += self.rho;             
        }
      }
    }
  } 


 fn observer(&self, x: Time, y: &State){// sample strains
      let mut sampled_strains = vec![0.0; self.nstrain];
      for k in 0..self.nhost {
         let y_k = y.slice_range(k*self.nstrain..(k+1)*self.nstrain,..);
         if y_k.sum() > 0.0 {
           let strain_dist = WeightedIndex::new(&y_k).unwrap();
           let s = strain_dist.sample(&mut rand::thread_rng());
           sampled_strains[s] += 1.0/(self.nhost as f64);
          }
       } 
        print!("{}", x);
        for s in sampled_strains {
          print!(",{}", s); 
        }
        print!("\n"); 
  }


/* ALTERNATE OBSERVER
  fn observer(&self, x: Time, y: &State){//return trajectories on host 0
      print!("{}", x);
      let k=0; // pick a host 
      let  y_k = y.slice_range(k*self.nstrain..(k+1)*self.nstrain,..);
      for y_i in &y_k {
        print!(",{}", *y_i); 
      }        
        print!("\n"); 
  }
*/
}


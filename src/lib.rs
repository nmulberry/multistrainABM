extern crate ode_event_solvers;
use ode_event_solvers::*;
extern crate rand;
extern crate rand_distr; 
extern crate array_tool;
use rand_distr::{Bernoulli, Poisson, Distribution};
use rand::distributions::WeightedIndex;
use rand::seq::index::sample;
use array_tool::vec::Intersect;
//========================================================//
/*                 MODEL DEFINITION                      */
//========================================================//
#[derive(Clone, Debug)]
pub struct WhCompetition {
  pub nhost: usize,
  pub nstrain: usize,
  pub theta: f64,
  pub eps_x: f64,
  pub eps_a: f64,
  pub rho: f64,
  pub alpha: Vec<f64>,
  pub kappa: Vec<f64>,
  pub id_m: Vec<u64>,
  pub id_s: Vec<u64>,
  pub id_o: Vec<u64>,
  pub beta: f64,
  pub t_g: f64,
  pub t_l: f64,
}

type State = DVector<f64>;
type Time = f64;

// Within Host Model for arbitrary number of strains and hosts
impl ode_event_solvers::System<State> for WhCompetition {
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

          // immunity
          dy[jj] = s[j]*s[j]/(self.rho*self.rho+s[j]*s[j])-
                     self.eps_a*y[jj]+
                     y[jj]*y[jj]/(self.theta*self.theta+
                                  y[jj]*y[jj]);
      }
    }
  }


  fn event(&mut self, _: Time, y: &State, dy: &mut State) {
    /* Here we model:
      - Transmission Events
      - Transformation Events (gain/loss)
      - Clearance Events (necessary fudge because of continous model)
      - Could add host level events here too
    */
    let trans_gain_dist = Bernoulli::new(self.t_g).unwrap();
    let trans_loss_dist = Bernoulli::new(self.t_l).unwrap();
    let caps_switch_dist = Bernoulli::new(self.t_g).unwrap();

    let mut n_infected:Vec<f64> = vec![0.0; self.nstrain];
    for k in 0..self.nhost {
      // get carriage
      for j in 0..self.nstrain {
        let ii = k*self.nstrain+j; //density of strain j
        if y[ii] < self.rho {
          dy[ii] = -y[ii];
        } else {        
           // track total density of the strain
          n_infected[j] += y[ii]; 
        }
      }
    }
    for j in 0..self.nstrain {
      if n_infected[j] > 0.0 {
        let transmit_dist = Poisson::new(self.beta*n_infected[j]).unwrap();
        let n_transmit = transmit_dist.sample(&mut rand::thread_rng());
        let recip = sample(&mut rand::thread_rng(), self.nhost, n_transmit as usize);
        for r in recip {
           dy[r*self.nstrain+j] += self.rho;             
        }
      }
    }

  // TRANSFORMATIONS
    for k in 0..self.nhost {
      // get carriage
      for j in 0..self.nstrain {
        let ii = k*self.nstrain+j; //density of strain j
        if y[ii] > self.rho {
          // first look for co-colonizations (have to loop over j again)
          for j2 in 0..self.nstrain {
             if j2 != j {
               let ii2 = k*self.nstrain+j2; //density of strain j2
               if y[ii2] > self.rho {
                  // co-colonization 
                  if self.id_s[j] != self.id_s[j2] {
                  // potential capsular switch
                  let caps_switch = caps_switch_dist.sample(&mut rand::thread_rng());
                  if caps_switch {
                    // create new recombinant strain in this host
                    // find which strain(note: wont be transmitted this day)
                    // sero from j, all other traits from j2
                    let ind1 = self.id_s
                          .iter()
                          .enumerate()
                          .filter(|(_, &r)| r == self.id_s[j])
                          .map(|(index, _)| index)
                          .collect::<Vec<_>>();
                    let ind2 = self.id_m
                          .iter()
                          .enumerate()
                          .filter(|(_, &r)| r == self.id_m[j2])
                          .map(|(index, _)| index)
                          .collect::<Vec<_>>();
                    let ind3 = self.id_o
                          .iter()
                          .enumerate()
                          .filter(|(_, &r)| r == self.id_o[j2])
                          .map(|(index, _)| index)
                          .collect::<Vec<_>>();

                    let j3 = ind3.intersect(ind1.intersect(ind2))[0];// SHOULD ONLY HAPPEN ONCE**
                    dy[k*self.nstrain+j3] += self.rho; // add new strain
                  }
                }
              }
            }
          }
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


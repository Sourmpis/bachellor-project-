from numpy import *
import nest
import nest.raster_plot
import pylab



def avg_firing_rate(spikes, dt, binsize, Tsim, Nneurons):
    """
      Calculates the average firing rate of a set of spike trains.
      spikes...all spike times of the population in units of [s]
      dt - A value of the firing rate is calculated at each time step dt [s]
      binsize - the size of the bin in multiples of dt. Used to calculate the
      firing rate at a particular moment.
                The rate is <num of spikes in [t-binsize*dt,t]> / (binsize*dt*Nneurons).
      Tsim - The length (in [s]) of the spike trains.
    """
    spi = array(floor(spikes/dt),dtype=int)
    Nbins = int(ceil(Tsim/dt))
    rate = zeros(Nbins+binsize-1)
    for sp in spi:
       rate[sp:sp+binsize]+=1
    rate /= (dt*binsize*Nneurons)
    rate = rate[:Nbins]
    return rate

def main():
    #general parameters
    Tsim = 1000.
    N_neurons = 1000
    N_E = int(0.8*N_neurons)
    N_I = int(0.2*N_neurons)
    C_E=int( N_E/2)
    C_I=int( N_I/2)


    #neuron parameters

    tau_m = 20.0
    E_L = -74.
    V_m =
    t_ref = 2.
    V_th = -54.0
    V_reset = -60.
    C_m = 10000.
    tau_syn_ex = tau_m
    tau_syn_in = tau_m
    E_ex = 0.
    E_in = -80.




    #Synapses parameters
    tau_psc = 3.
    tau_fac = 41.
    tau_rec = 26.
    U = 0.5
    delay_syn = 0.1
    u = 0.
    x = 1.
    g = 5.
    J_E   = 1000.
    g = 5.0
    J_I   =  -g * J_E

    #spike generator parameters
    eta = 20.0
    nu_ex   = eta * (V_th-E_L) /( J_E * C_E * tau_m )
    p_rate = 1000.0 * nu_ex  * C_E

    syn_param = {"tau_psc": tau_psc,
             "tau_fac": tau_fac,  # facilitation time constant in ms
             "tau_rec": tau_rec,  # recovery time constant in ms
             "U": U,             # utilization
             "delay": delay_syn,       # transmission delay
             "weight": J_E,
             "u": u,
             "x": x}


    nrn_params =     {    "V_m":V_rest_ex,     # Membrane potential in mV
                          "E_L": E_L,     # Resting membrane potential in mV
                          "C_m": C_m,           # Capacity of the membrane in pF
                          "tau_m": tau_m,      # Membrane time constant in ms, R_m*tau_m
                          "V_th": V_th,     # Spike threshold in mV
                          "V_reset": V_reset, # Reset potential of the membrane in mV
                  "t_ref": 1.   # refractory time in ms
                          }


    nest.ResetKernel()
    nest.SetKernelStatus({"resolution": 0.1})

    nodes_E = nest.Create('iaf_psc_exp',N_E,nrn_params_ex)
    nodes_I = nest.Create('iaf_psc_exp',N_I,nrn_params_in)
    nodes = nodes_E+nodes_I
    nest.CopyModel("tsodyks_synapse",
                   'excitatory',
                   syn_param)

    nest.Connect(nodes_E,nodes,
                 {'rule': 'fixed_indegree', 'indegree': C_E},
                 syn_spec = "excitatory")

    nest.CopyModel('tsodyks_synapse',
                   'inhibitory',
                   syn_param)

    nest.Connect(nodes_I,nodes,
                 {'rule': 'fixed_indegree','indegree': C_I},
                 syn_spec = 'inhibitory')

    noise = nest.Create('poisson_generator',20,{'rate':p_rate,'start':0.0001})

    nest.Connect(noise,nodes_E[:20],{'rule':'one_to_one'},syn_spec={'weight':10000000.})

    spikes = nest.Create('spike_detector',2,[{'label':'brunel-py-ex'},
                                             {'label': 'brunel-py-in'}])
    spikes_E = [spikes[0]]
    spikes_I = [spikes[1]]
    N_rec = 100
    nest.Connect(nodes_E,spikes_E)
    nest.Connect(nodes_I,spikes_I)

    nest.Simulate(200.)

    print('fuck pynest')

    nest.GetStatus(noise)
    nest.SetStatus(noise,{'rate':0.})

    nest.Simulate(Tsim-200.)


    events = nest.GetStatus(spikes,keys = "events")
    evs_ex = events[0]["senders"]
    ts_ex = events[0]["times"]
    evs_in = events[1]["senders"]
    ts_in = events[1]["times"]


    rate_ex = len(ts_ex)/Tsim/N_rec
    print("Excitatory rate :" ,rate_ex)
    rate_in = len(ts_in)/Tsim*1000./N_rec
    print("Inhibitory rate :" ,rate_in)

    dt = 0.0005
    binsize = 10

    rate_ex_1 = avg_firing_rate(ts_ex/1000., dt, binsize, Tsim/1000., N_neurons)


    pylab.figure(1)


    t = linspace(0,0.1,2000)
    print(shape(t),shape(rate_ex_1))
    pylab.plot(t,rate_ex_1)

    # nest.raster_plot.from_device([spikes[0]],hist = True)
    # nest.raster_plot.from_device([spikes[1]],hist = True)

    pylab.show()

if __name__ == "__main__":
    main()

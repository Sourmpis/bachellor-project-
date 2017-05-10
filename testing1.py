import nest
import nest.raster_plot
import pylab
Tsim = 300.
g = 5.0
eta = 20.0
delay = 1.5
tau_m = 20.0
C_m = 1.
t_ref = 2.
V_th = 20.0
V_reset = 10.0
V_rest = 0.
N_E = 8000
N_I = 2000
N_neurons = N_E+N_I
C_E=int( N_E/10)
C_I=int( N_I/10)
J_E   = 1000.
J_I   =  -g * J_E
nu_ex   = eta * V_th /( J_E * C_E * tau_m )
p_rate = 1000.0 * nu_ex  * C_E

nest.ResetKernel()
nest.SetKernelStatus({"resolution": 0.1})
nest.SetDefaults('iaf_psc_delta',
                 {'C_m': C_m,
                  'tau_m': tau_m,
                  't_ref' : t_ref,
                  'V_th' : V_th,
                  'V_reset': V_reset,
                  'V_m':V_rest
                 })
nodes = nest.Create('iaf_psc_delta',N_neurons)
nodes_E = nodes[:N_E]
nodes_I = nodes[N_E:]

nest.CopyModel('static_synapse_hom_w',
               'excitatory',
               {'weight':J_E,
                'delay':delay})

nest.Connect(nodes_E,nodes,
             {'rule': 'fixed_indegree', 'indegree': C_E},
             syn_spec = "excitatory")

nest.CopyModel('static_synapse_hom_w',
               'inhibitory',
               {'weight':J_I,
                'delay':delay})

nest.Connect(nodes_I,nodes,
             {'rule': 'fixed_indegree','indegree': C_I},
             syn_spec = 'inhibitory')

noise = nest.Create('poisson_generator',1,{'rate':p_rate})

nest.Connect(noise,nodes,syn_spec = 'excitatory')

spikes = nest.Create('spike_detector',2,[{'label':'brunel-py-ex'},
                                         {'label': 'brunel-py-in'}])
spikes_E = [spikes[0]]
spikes_I = [spikes[1]]
N_rec = 50
nest.Connect(nodes_E[:N_rec],spikes_E)
nest.Connect(nodes_I[:N_rec],spikes_I)

nest.Simulate(Tsim/2)

print('fuck pynest')

nest.GetStatus(noise)
nest.SetStatus(noise,{'rate':0.})

nest.Simulate(Tsim/2)


events = nest.GetStatus(spikes,keys = "events")
evs_ex = events[0]["senders"]
ts_ex = events[0]["times"]
evs_in = events[1]["senders"]
ts_in = events[1]["times"]


rate_ex = len(ts_ex)/Tsim*1000./N_rec
print("Excitatory rate :" ,rate_ex)
rate_in = len(ts_in)/Tsim*1000./N_rec
print("Inhibitory rate :" ,rate_in)



nest.raster_plot.from_device([spikes[0]],hist = True)
nest.raster_plot.from_device([spikes[1]],hist = True)

pylab.show()
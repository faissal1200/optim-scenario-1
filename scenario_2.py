from gurobipy import *
import numpy as np
import json
import os
import matplotlib.pyplot as plt
from matplotlib import gridspec

end = 35121
#end = 96 * 7
data_folder = 'profils'
cooling_file = 'aggregated_cooling_needs.json'
with open(os.path.join(data_folder, cooling_file), 'r') as file:
    cooling_needs = json.load(file)[:end]

heat_needs_gas = 'aggregated_heating_needs_gas.json'
with open(os.path.join(data_folder, heat_needs_gas), 'r') as file:
    heating_needs_gas = json.load(file)[:end]

heat_needs_net = 'aggregated_heating_needs_net.json'
with open(os.path.join(data_folder, heat_needs_net), 'r') as file:
    heating_needs_net = json.load(file)[:end]

pv_prod = 'aggregated_pv_prod.json'
with open(os.path.join(data_folder, pv_prod), 'r') as file:
    solar_prod = json.load(file)[:end]

sub_heat = 'aggregated_heating_needs_sub.json'
with open(os.path.join(data_folder, sub_heat), 'r') as file:
    sub_need = json.load(file)[:end]

elec_direct = "aggregated_elec_services.json"
with open(os.path.join(data_folder, elec_direct), 'r') as file:
    electricity_need = json.load(file)[:end]


np.random.seed(0)

m = Model('energy_system')

r_p2g = 0.55  # low Voltage Power Grid to power to gas
r_gas_b = 0.9  # ratio gas_network to gas_boiler
r_chp_lv = 0.3  # ratio combined to low Voltage Power Grid
r_chp_ht = 0.6  # ratio combined heat pump to high tension network
r_co2_gas = 0.184  # co2 emissions/kwh of gas
r_co2_pv = 0.09  # co2 emissions/kwh of solar energy
r_co2_mv = 0.0236  # co2 emissions/kwh of solar energy
r_co2_chp = 0.0236  # co2 emissions/kwh of solar energy


# Creating needs profile
def spread(seq, n=4):
    seq_res = np.array([])
    for i in range(1, len(seq)):
        spr = np.linspace(seq[i - 1], seq[i], num=n, endpoint=False)
        seq_res = np.concatenate((seq_res, spr))
    return seq_res


t_heat = 60  # heating heat pump temperature[deg.C]
t_cool = 2  # cooling heat pump temperature[deg.C]
t_lake = 8  # lake temperature[deg.C]
t_res_1 = 25  # LT network temperature [deg.C]
t_heat_2 = 80  # mid heat pump temperature[deg.C]


# COP for heating heat pump computing
def cop_heat(t_con, t_eva):
    return 0.4 * (t_con + 273.15) / (t_con - t_eva)


# COP for cooling heating pump computing
def cop_cool(t_con, t_eva):
    return 0.4 * (t_eva + 273.15) / (t_con - t_eva)


cop = {('hp_lake', t_res_1): cop_heat(t_res_1, t_lake),
       ('hp_heat_lt', t_res_1): cop_heat(t_heat, t_res_1),
       ('hp_cool_lt', t_res_1): cop_cool(t_res_1, t_cool), ('hp_heat_ht', t_res_1): cop_heat(t_heat_2, t_res_1)}

base_links = [
    # lt_dhcn
    ('hp_lake', 'lt_dhcn'),
    ('lake', 'hp_lake'),
    ('hp_cool_lt', 'lt_dhcn'),
    ('lt_dhcn', 'hp_heat_lt'),
    ('lt_dhcn', 'reject'),
    ('lt_dhcn', 'hp_heat_ht'),

    # ht_dhcn
    ('hp_heat_ht', 'ht_dhn'),
    ('ht_dhn', 'sub'),
    ('gas_b', 'ht_dhn'),
    ('chp', 'ht_dhn'),

    # Gas
    ('p2g', 'gas'),
    ('gas', 'chp'),
    ('gas', 'gaz_direct'),
    ('gas', 'gas_b'),
    ('mp', 'gas'),

    # input_LV
    ('pv', 'lv'),
    ('chp', 'lv'),
    ('mv', 'lv'),

    # output_LV
    ('lv', 'hp_heat_ht'),
    ('lv', 'hp_cool_lt'),
    ('lv', 'hp_heat_lt'),
    ('lv', 'hp_lake'),
    ('lv', 'electricity_direct'),
    ('lv', 'p2g'),
    ('lv', 'mv'),

]

base_links_2 = [
    ('gas_storage', 'gas'),

]

links_1 = tuplelist([(l[0], l[1]) for l in base_links])
links_2 = tuplelist([(l[0], l[1]) for l in base_links_2])

var_links = m.addVars(links_1, end, lb=0)  # Creating link for each time_step for base_links
var_links_2 = m.addVars(links_2, end, lb=-GRB.INFINITY)  # Creating link for each time_step for base_links

# Creation of the cumulative variable
cum_sum = m.addVars(end)
mv_var_max = var_links.select('lv', 'mv', '*')
mv_var_min = var_links.select('lv', 'mv', '*')
normal = var_links_2.select('gas_storage', 'gas', '*')

m.addConstrs(cum_sum[i] == normal[i] + cum_sum[i - 1]
             for i in range(1, end))
m.addConstr(cum_sum[0] == normal[0])

value_max = m.addVar()
value_min = m.addVar()
value_mv_max = m.addVar()
value_mv_min = m.addVar()

m.addGenConstrMin(value_min, cum_sum)  # minimum of the cumulative variable
m.addGenConstrMax(value_max, cum_sum)  # maximum of the cumulative variable
m.addGenConstrMax(value_mv_max, mv_var_max)  # maximum of the cumulative variable
m.addGenConstrMin(value_mv_min, mv_var_min)  # minimum of the cumulative variable

var_links.update(var_links_2)


# Creation of conservation constraint for each hub
def hub(net):
    m.addConstrs(
        (var_links.sum('*', net, idt) == var_links.sum(net, '*', idt)
         for idt in range(end)))


hub('ht_dhn')
hub('lt_dhcn')
hub('lv')
hub('gas')

# gas storage
m.addConstr(
    var_links.sum('gas_storage', 'gas', '*') == 0,
    name='gas_storage')

# ht, sub_station
m.addConstrs(
    (var_links.sum('ht_dhn', 'sub', idt) == sub_need[idt] / 1000
     for idt in range(end)))

# Solar production
m.addConstrs(
    (var_links.sum('pv', 'lv', idt) == solar_prod[idt] * (100 / 1000)
     for idt in range(end)))


# gas, gas_direct
m.addConstrs(
    (var_links.sum('gas', 'gaz_direct', idt) == heating_needs_gas[idt] / 1000
     for idt in range(end)))


# lv, electricity_demand
m.addConstrs(
    (var_links.sum('lv', 'electricity_direct', idt) == electricity_need[idt] / 1000
     for idt in range(end)), name='cop heating')

# hp_cool_lt, lt_dhcn
m.addConstrs(
    (var_links.sum('hp_cool_lt', 'lt_dhcn', idt) == (1 + 1 / (cop['hp_cool_lt', t_res_1])) * (cooling_needs[idt] / 1000)
     for idt in range(end)), name='cop cooling')

# low voltage grid, hp_cool_lt
m.addConstrs(
    (var_links.sum('lv', 'hp_cool_lt', idt) == (cooling_needs[idt] / 1000) / cop['hp_cool_lt', t_res_1]
     for idt in range(end)), name='cop cooling')

# lt_dhcn, hp_heat_lt
m.addConstrs(
    (var_links.sum('lt_dhcn', 'hp_heat_lt', idt) == (1 - 1 / (cop['hp_heat_lt', t_res_1])) * (heating_needs_net[idt] / 1000)
     for idt in range(end)), name='cop heating')

# low voltage grid, hp_heat_lt
m.addConstrs(
    (var_links.sum('lv', 'hp_heat_lt', idt) == (heating_needs_net[idt] / 1000) / cop['hp_heat_lt', t_res_1]
     for idt in range(end)), name='cop heating')

# hp_heat_ht, ht_dhn
m.addConstrs(
    (var_links.sum('hp_heat_ht', 'ht_dhn', idt) == (1 + 1 / cop['hp_heat_ht', t_res_1]) * var_links.sum('lt_dhcn',
                                                                                                        'hp_heat_ht',
                                                                                                        idt)
     for idt in range(end)), name='cop heating')
m.addConstrs(
    (var_links.sum('lv', 'hp_heat_ht', idt) * cop['hp_heat_ht', t_res_1] == var_links.sum('hp_heat_ht', 'ht_dhn', idt)
     for idt in range(end)), name='cop heating')

# low voltage grid, hp_lake
m.addConstrs((var_links.sum('lv', 'hp_lake', idt) * cop['hp_lake', t_res_1] == var_links.sum('hp_lake', 'lt_dhcn', idt)
              for idt in range(end)), name='cop electric lake')

# p2g, gas
m.addConstrs(
    (var_links.sum('p2g', 'gas', idt) == r_p2g * var_links.sum('lv', 'p2g', idt)
     for idt in range(end)))


# gas, gas_b, ht_dhn
m.addConstrs(
    (var_links.sum('gas_b', 'ht_dhn', idt) == r_gas_b * var_links.sum('gas', 'gas_b', idt)
     for idt in range(end)))
# chp, ht_dhn
m.addConstrs(
    (var_links.sum('chp', 'ht_dhn', idt) == r_chp_ht * var_links.sum('gas', 'chp', idt)
     for idt in range(end)))
# chp, lv
m.addConstrs(
    (var_links.sum('chp', 'lv', idt) == r_chp_lv * var_links.sum('gas', 'chp', idt)
     for idt in range(end)))

# lake, hp_lake
m.addConstrs(
    (var_links.sum('lake', 'hp_lake', idt) == var_links.sum('hp_lake', 'lt_dhcn', idt) - var_links.sum('lv', 'hp_lake',idt)
     for idt in range(end)))

# creation of the objective function elements
feed_electricity_mv = var_links.sum('mv', 'lv', '*')
feed_gas_mp = var_links.sum('mp', 'gas', '*')
ext_feeder = m.addVar(name='ext_feeder')

m.addConstr(ext_feeder == 17.68 * feed_electricity_mv + 9.5 * feed_gas_mp, name='ext_feeder')

# multi objective functions
m.ModelSense = GRB.MINIMIZE
m.NumObj = 2

m.setParam(GRB.Param.ObjNumber, 0)
m.ObjNPriority = 2
m.ObjNWeight = 1
m.ObjNName = 'Ext_feeder'
m.ObjNRelTol = 0.1
ext_feeder.ObjN = 1.0

m.setParam(GRB.Param.ObjNumber, 1)
m.ObjNPriority = 1
m.ObjNWeight = 1.0
m.ObjNName = 'storage'
value_max.ObjN = 1.0
value_min.ObjN = -1.0

m.optimize()

pv = np.array([v.x for v in var_links.select('pv', 'lv', '*')])
mv = np.array([v.x for v in var_links.select('mv', 'lv', '*')])
lv_mv = np.array([v.x for v in var_links.select('lv', 'mv', '*')])

stck = np.array([v.x for v in var_links.select('gas_storage', 'gas', '*')])
fig = plt.figure(0)
gs = gridspec.GridSpec(2, 2)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])

ax1.plot(pv, '.', color='green', alpha=0.5, label='pv')
ax1.plot(mv, '.', color='blue', alpha=0.5, label='mv')
ax1.plot(lv_mv, '.', color='red', alpha=0.5, label='lv_mv')
ax2.plot(stck, '.', color='blue', alpha=0.5, label='storage')
ax1.legend()
ax2.legend()

plt.figure(1)
plt.plot(mv, '.', color='blue', alpha=0.5, label='MV vers LV')
plt.plot(lv_mv, '.', color='red', alpha=0.5, label='LV vers MV')
plt.ylabel(' Electricité [kw]')
plt.xlabel('Temps en 1/4 h')
plt.legend()

plt.figure(2)
plt.plot(pv, '.', color='green', alpha=0.5, label='PV')
plt.plot(mv, '.', color='blue', alpha=0.5, label='MV')
plt.ylabel('Electricité [kw]')
plt.xlabel('Temps en 1/4 h')
plt.legend()

plt.figure(3)
plt.plot(stck, '.', color='blue', alpha=0.5, label='stockage')
plt.ylabel('gaz')
plt.xlabel('Temps en 1/4 h')
plt.legend()
plt.show()

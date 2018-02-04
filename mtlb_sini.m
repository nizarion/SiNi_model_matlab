clear;
%%
%parameters
%constants
%[a.u.]       [J/K]       [C]	         [F/m]
pi=3.1415926; Kb=1.38e-23; q=1.6e-19; eps_0=8.85e-12;
%conductance parameters
%[m^2]      [m]     [eV]
area=1e-8; d=7e-9; Ub=0.85;
%+ area=6.4e-9 d=4e-8  Ub=0.9
%[1/V]    [1/V]  [A/(m^2 K^2)] [A/(m^2 K^2)]
scl_p=5;  scl_n=7; AA_p=2.4e6;   AA_n=1e6;
%+ scl_p=5  scl_n=7 AA_p=2.4e6   AA_n=1e6
% [a.u.] [ohm]   [K]
eps_r=5.5; Ron=300; T0=300; k_temp=1;
% + eps_r=8 Ron=300 T0=300 k_temp=1
% [A]  [ohms]      [ohms]
cI =1e-3; r_on_p=20e3; r_on_n=30e3;
%+ cI=1e-4 r_on_p=6.1k r_on_n=8.5k
% [F/m]
eps_i=eps_r*eps_0;
%energy threshold parameters
p_th_set=30e-11; p_th_reset=-4e-5;
% p_th_set=1e-10 p_th_reset=-3e-7
% windowing: soft switching
% + p=4
% voltage dependent variabiligy: windowing
p_rand_i=4;
%+ p_rand_i=4
% internal params
v_off=0; v_on=1;
%+ v_off=0 v_on=1e-7
%Poole-Frenkel constant
pfc=1e-1;
%%
%Initialazation

h=0.001;	%the step of runge-kutta method (time precision)  
k=1; %index initialization
t=0; %time initialization
time(k)=0;
Im(k)=0;
Vm(k)=0;		 %initialization of the memristor's voltage

%%
%Programming/bias signal amplitude and fequency.
Vo=5;		%amplitude of the input voltage
freq=1;     %frequency of the generated pulses 
%T=1/freq;
%%
%Initial Conditions
p_sum_pow=0;
n_sum_pow=-40.139e-6;
set_th=p_th_set;
reset_th=p_th_reset;
state(k)=v_off;

i_G_off_p = 0;
i_G_off_n = 0;
%%
while (t<=2)
    %%
    %Programming signal declaration:
    Vm(k)=Vo*sin(2*pi*freq*t);
    %%
    %energy computation
    %v= idt( v(Plus, Minus)>0 && v(p_sum_pow)<v(set_th) ? v(Plus, Minus*abs(i(G_cond)) : 0, 0 )
    
    if (Vm(k) > 0) && (p_sum_pow<set_th)
        p_sum_pow = p_sum_pow + ( Vm(k)*abs(Im(k)) * h );
    else
        %p_sum_pow = p_sum_pow + 0 ;
    end
    %n_sum_pow v= idt( v(Plus, Minus)<0 && v(n_sum_pow)>v(reset_th) ?
    %    +v(Plus, Minus)*abs(i(G_cond)) : 0, 0 )
    if (Vm(k)< 0) && (n_sum_pow>reset_th)
        n_sum_pow = n_sum_pow + ( Vm(k)*abs(Im(k)) * h );
    else
        %p_sum_pow = p_sum_pow + 0 ;
    end
    
    %%
    %thresholds computation
    %v= ( floor( v(n_sum_pow)/p_th_reset ) + 1 )*p_th_set
    
    set_th = (floor( n_sum_pow / p_th_reset) +1) * p_th_set;
    reset_th = (floor( p_sum_pow / p_th_set) +1) * p_th_reset;
    
    %%
    %state computation
    if state(k) == v_off
        if p_sum_pow>set_th
            state(k+1) = v_on;
        else
            state(k+1) = v_off;
        end
    elseif state(k) == v_on
        if n_sum_pow<reset_th
            state(k+1) = v_off;
        else
            state(k+1) = v_on;
        end
    else
        errror = 1 ;
    end
        %%
        %conduction computation
        E_on = Vm(k);

        i_Ron_p = E_on/ r_on_p;
        i_Ron_n = E_on/ r_on_n;

        %%
        %off modeling
        i_G_off_p = area* pfc * abs(Vm(k)) /d  *exp(-q * Ub / (Kb * T0) ) ...
            * exp( sqrt(abs(Vm(k)) ) * (q / (Kb * T0) * ...
            ( sqrt(q / (d * 4 * pi * eps_i)))) );

        i_G_off_n = -area* pfc * abs(Vm(k))/d  *exp(-q * Ub / (Kb * T0) ) ...
            * exp( sqrt(abs(Vm(k)) ) * (q / (Kb * T0) * ...
            ( sqrt(q / (d * 4 * pi * eps_i)))) );
    %%
    %Current source selection
    if state(k) == v_off
        if (Vm(k) > 0)
            Im(k+1) = i_G_off_p;
%             Vm(k+1)=i_G_off_p*r_on_p;
        else
            Im(k+1) = i_G_off_n;
%             Vm(k+1)=i_G_off_n*r_on_n;
        end
    else
        if (Vm(k) > 0)
            Im(k+1) = i_Ron_p;
%             Vm(k+1)=i_Ron_p*r_on_p;
        else
            Im(k+1) = i_Ron_n;
%             Vm(k+1)=i_Ron_n*r_on_n;
        end
        
    end
    Im(k+1) = min( cI, Im(k+1));
    Im(k+1) = max(-cI, Im(k+1));
    %%
    %cycles computation
    % variable cycle counts how many time the meristor's 
    %been throughon-off cycle
	cycle = floor( n_sum_pow/p_th_reset );
    %variable sEvent counts how many times an event e.g. switching to 
    %on state, sw to off state has occured
    sEvent = floor(  n_sum_pow/p_th_reset ) + floor(p_sum_pow/p_th_set );
    %%
    time(k)=t; %usefull for plotting the device responce (time vector)
    t = t+h; %time precision defined by the solver's step
    k = k+1; %increase index

end
Im=Im(1:end-1);
state=state(1:end-1);
%%
%Plots:

figure('Position', [10, 10, 1000, 750])
subplot(2,3,[2 3 5 6])
plot(Vm,log(abs(Im))),grid on, title('Pinched hysterisis loop (I-V characteristic)'),xlabel('voltage(Vm)'),ylabel('current(Im)');

% subplot(2,2,3)
% plot(time,L),grid on, title('Tunnel barrier width (L)'),xlabel('time(t)'),ylabel('barrier width(L)');
% 
% subplot(2,2,4)
% plot(time,RL),grid on, title('Memristance'),xlabel('time(t)'),ylabel('memristance(RL)');
subplot(2,3,1)
plot(time,state)
hold on;
subplot(2,3,1)
plot(time,Vm),grid on, title('Programming voltage (Vm)'),xlabel('time(t)'),ylabel('voltage(Vm)');
hold off;
subplot(2,3,4)
plot(time,Im),grid on,  title('Memristor current (Im)'),xlabel('time(t)'),ylabel('current(Im)');
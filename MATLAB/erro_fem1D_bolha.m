function [L2error,Eerror,Herror] = erro_fem1D_bolha(x,u,ufun,dufun)
% Devolve as respetivas normas
    % Valores para as normas
    val_L2=0;
    val_E=0;
    val_H=0;
    for i=1:(length(x)-1)
        % Regra de Quadratura de Gauss
        [t,w] = gaussint(3,x(i),x(i+1));
        % Pontos avaliados
        uh_tk = zeros(1,length(t));
        duh_tk = zeros(1,length(t));
        u_tk = zeros(1,length(t));
        du_tk = zeros(1,length(t));
            for k=1:length(t)
                uh_tk(k) = u(i)*( x(i+1)-t(k))/(x(i+1)-x(i));
                uh_tk(k) = uh_tk(k) + u(i+1)*( x(i)-t(k))/(x(i)-x(i+1));
                uh_tk(k) = uh_tk(k) + u(i+length(x))*((t(k)-x(i))*(x(i+1)-t(k)));
                duh_tk(k) = -u(i)/(x(i+1)-x(i));
                duh_tk(k) = duh_tk(k) - u(i+1)/(x(i)-x(i+1));
                duh_tk(k) = duh_tk(k) + u(i+length(x))*((-2*t(k))+x(i)+x(i+1));
                u_tk(k) = ufun(t(k));
                du_tk(k) = dufun(t(k));
            end
        val_L2 = val_L2 + (u_tk - uh_tk).^2*w(:);
        val_E = val_E + (du_tk - duh_tk).^2*w(:);
        val_H = val_H + (u_tk - uh_tk).^2*w(:) + (du_tk - duh_tk).^2*w(:) ;
    end
    L2error = sqrt(val_L2);
    Eerror = sqrt(val_E);
    Herror = sqrt(val_H);
end
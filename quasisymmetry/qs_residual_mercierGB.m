function [f,Jacf] = qs_residual_mercierGB(sigma,mu,delta,etab,curv)
    
    f = [sigma-mu.*sin(2.*delta)./sqrt(1-mu.^2);
         etab^2-((1-mu.*cos(2.*delta))./sqrt(1-mu.^2)).*(curv.^2)];
     
    df1dmu    = -diag(((mu.^2).*sin(2.*delta))./((1-mu.^2).^(1.5))+sin(2.*delta)./sqrt(1-mu.^2));
    df1ddelta = diag(-2.*mu.*cos(2.*delta)./sqrt(1-mu.^2));
    df2dmu    = diag(cos(2.*delta)./sqrt(1 - mu.^2) - (mu.*(1 - mu.*cos(2.*delta)))./((1 - mu.^2).^(1.5)));
    df2ddelta = diag(-2.*mu.*sin(2.*delta)./sqrt(1-mu.^2));

    Jacf = [df1dmu df1ddelta;
            df2dmu df2ddelta];

end
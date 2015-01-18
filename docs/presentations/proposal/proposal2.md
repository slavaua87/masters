---
title       : Adaptive Action Selection
subtitle    : Master's Proposal
author      : Slava Nikitin
framework   : io2012
highlighter : highlight.js  # {highlight.js, prettify, highlight}
hitheme     : tomorrow    
widgets     : [mathjax] 
mode        : standalone
knit        : slidify::knit2slides

--- 
<!-- Limit image width and height -->
<style type="text/css">
img {    
  max-height: 560px;    
  max-width: 964px;
}
</style>
 
<!-- Center image on slide -->
<script type="text/javascript" src="http://ajax.aspnetcdn.com/ajax/jQuery/jquery-1.7.min.js"></script>
<script type="text/javascript">
$(function() {    
  $("p:has(img)").addClass('centered');
});
</script>

## Psychological Problem

* Action selection guided by noisy, analogue sensory signal

--- &twocol wd1:50% wd2:50%
## Ratcliff & Rouder (1998) Task

*** =left
![Bright](https://dl.dropboxusercontent.com/u/14983136/bright_patch.jpeg)

*** =right
![Dark](https://dl.dropboxusercontent.com/u/14983136/dark_patch.jpeg)


--- 
## Theoretical Targets

* Decision mechanism and its across-trial adaptation
  + Thresholded accumulators
  + Stimulus quality and probability, speed-accuracy settings, confidence

--- 
## Decision Model

With drift coefficient $\delta \in \mathbb{R}$ and diffusion coefficient $\sigma > 0$, Wiener process starts at $x(0) = \xi$ and evolves over state space $\mathcal{S} = [0, \alpha]$ until it is absorbed at one of the boundaries as follows:

$$
dx(t) = \delta dt + \sigma\sqrt{dt}dw(t)
$$

![Accumulator](https://dl.dropboxusercontent.com/u/14983136/diffusion.jpeg)

--- 
## Modeling Adaptation

* Adaptation as variation of parameters across trials.
  + Dynamical approach:
$$
\begin{eqnarray}
\xi(n+1) = f_{\xi}(\xi(n), t^{rt}(n), r(n),\theta_{\xi}) \\
\delta(n+1) = f_{\delta}(\delta(n), t^{rt}(n), r(n),\theta_{\delta}) \\
\end{eqnarray}
$$
  + Statistical approach:
$$
\begin{equation}
(\xi(n+1), \delta(n+1)) \sim f(\theta_{\xi}, \theta_{\delta} \mid \xi(n), \delta(n))
\end{equation}
$$

--- 
## Statistical Approach

* Parameter variation structure:
  + Auto-correlation
  + Cross-correlation

---
## Modeling Cross-correlation 

* $\textbf{Definition}$:
A function $C(u_1, u_2, \ldots, u_p \mid \boldsymbol{\theta}): [0, 1]^p \mapsto [0, 1]$ is a copula with parameter $\boldsymbol{\theta} \in \Theta \subseteq \mathbb{R}^k$ if $C$ is a multivariate probability distribution function with univariate marginals $u_i$ distributed as standard uniform random variables.  

---
## Copula Theory

* $\textbf{Sklar's theorem}$:  
Let $F$ be a multivariate distribution with continuous marginal
distributions $F_1, F_2, \dots, F_p$, then there exists a unique
copula $C$ such that 
$$
\begin{equation}
F(x_1, x_2, \dots, x_p \mid \boldsymbol{\psi}, \boldsymbol{\theta}) = C(F_1(x_1 \mid \psi_1), F_2(x_2 \mid \psi_2), \dots, F_p(x_p \mid \psi_p) \mid \boldsymbol{\theta}).
\end{equation}
$$
Conversely, given a copula $C$ and marginal univariate distributions $F_1, F_2, \dots, F_p$, $F$ is a multivariate distribution.

* $\textbf{Corollary}$:  
$$
\begin{equation}
f(x_1, x_2, \dots, x_p \mid \boldsymbol{\psi}, \boldsymbol{\theta}) = \\ 
c(F_1(x_1 \mid \psi_1), F_2(x_2 \mid \psi_2), \dots, F_p(x_p \mid \psi_p) \mid \boldsymbol{\theta})\prod_{i = 1}^p f(x_i \mid \psi_i)
\end{equation}
$$

---
## Standard Diffusion Model Example

Assume normal distribution $F_n$ for the drift $\delta$ and uniform distribution $F_u$ for the starting point $\xi$, then
$$
\begin{equation}
F(\delta, \xi) = F_n(\delta)F_u(\xi),
\end{equation}
$$
and the copula for the standard diffusion model is 
$$
\begin{equation}
C(u_1,u_2) = u_1 u_2.
\end{equation}
$$

![diff_stand](https://dl.dropboxusercontent.com/u/14983136/diff_stand.jpeg)

---
## Example with Correlated Parameters

For positive correlation, use a Frank copula 
$$
\begin{equation}
C(u_1, u_2) = -\frac{1}{\theta}\ln \left(1 + \frac{\left(\exp(-\theta u_1)-1\right)\left(\exp(-\theta u_2)-1\right)}{\left(\exp(-\theta)-1\right)}
\right),
\end{equation}
$$
where $\theta_F \in \mathbb{R} \setminus \{0\}$ is a dependency parameter. Then a bivariate distribution of $(\delta, \xi)$ is
$$
\begin{equation}
F(\delta, \xi \mid \theta, \mu, \sigma^2, a, b) = C(F_n(u_1 \mid \mu, \sigma^2), F_u(u_2 \mid a, b) \mid \theta).
\end{equation}
$$

![diff_cor](https://dl.dropboxusercontent.com/u/14983136/diff_copulas.jpeg)

---
## Parameter Model Desiderata

1. Copulas have closed form
2. Each dimension has a separate correlation parameter
3. Correlation spans from -1 to 1, or full range of dependence
4. Closed form marginals with flexible shape

--- 
## Copulas
$\textbf{Independence}$
$$
\begin{equation}
C(u_1, u_2, u_3) = u_1 u_2 u_3
\end{equation}
$$
$\textbf{Normal}$
$$
\begin{equation}
C(u_1, u_2, u_3 \mid \boldsymbol{P}_n) = \int_{-\infty}^{\Phi^{-1}(u_1)} \int_{-\infty}^{\Phi^{-1}(u_2)} 
        \int_{-\infty}^{\Phi^{-1}(u_3)} 
        |(2\pi\boldsymbol{P}_n)|^{-1 / 2}\exp\left(-\frac{1}{2}\boldsymbol{x}^T
\boldsymbol{P}_n \bf{x}\right)d\bf{x}
\end{equation}
$$
$\textbf{t}$
$$
\begin{equation}
C(u_1, u_2, u_3 \mid \omega, \boldsymbol{P}_t) = \int_{-\infty}^{t^{-1}(u_1 \mid \omega)}
       \int_{-\infty}^{t^{-1}(u_2 \mid \omega)}
       \int_{-\infty}^{t^{-1}(u_3 \mid \omega)}
       \left(1 + \boldsymbol{x}^T\boldsymbol{P}_t^{-1 / 2}\boldsymbol{x} /
{\omega} \right)^{(\omega + 3) / 2} \\
       \times \mathcal{G} \left([\omega+3] / 2\right) / \left[\mathcal{G}(2\omega)\sqrt{|\pi\omega\boldsymbol{P}_t|}\right] d\boldsymbol{x}
\end{equation}
$$

---
## Univariate Marginals

Parameterize with bias $\beta = \xi / \alpha \in (0, 1)$ , so univariate marginals are
$$
\begin{eqnarray}
\delta & \sim & \mathcal{N}(\nu, \eta^2), \nonumber \\
\beta & \sim & \mathcal{B}(\lambda, \gamma), \text{and} \nonumber \\
t^{nd} & \sim & \mathcal{TN}(\chi,\phi^2,a,b).
\end{eqnarray}
$$
Combining copula $i \in \{1, 2\}$ and univariate marginals leads to a multivariate distribution
$$
\begin{equation}
F_i(\delta,\beta,t^{nd}) = C_i(F_1(\delta),F_2(\beta),F_3(t^{nd})),
\end{equation}
$$

or to its related density
$$
\begin{equation}
f_i(\delta,\beta,t^{nd}) 
   = c_i(F_1(\delta),
         F_2(\beta),
         F_3(t^{nd}))
         f_1(\delta)f_2(\beta)f_3(t^{nd}).
\end{equation}
$$

--- 
## Thesis Problems

* Theoretical

* Empirical

---
## Features of the Theoretical Study

1. 
$$
\begin{equation}
\mathrm{E}[x(t)] = \int_{\mathcal{S}} x(t) f(x(t)) dx(t)
\end{equation}
$$
2.  
$$
\begin{equation}
f(t^{rt} \mid r, \boldsymbol{\theta}) = \int_A f(t^{rt}, r \mid \alpha, \delta, \xi) f(\delta, \xi \mid \boldsymbol{\theta}) dA
\end{equation}
$$
3. Ratcliff and Rouder (1998) design with 9 combinations of correlations.

---
## Features of the Empirical Study

1. Outlier processing
2. Hierarchical Bayesian estimation 
3. Differential evolution MCMC
$$
\boldsymbol{\rho}^* = \boldsymbol{\rho}^{l, j-1} + \upsilon(\boldsymbol{\rho}^{m, j - 1} - \boldsymbol{\rho}^{n, j - 1}) + \boldsymbol{\varepsilon}
$$

---
## Analyses

1. Probabilities, HPD intervals
2. Posterior prediction checks
3. $BPIC = \operatorname{E_{\boldsymbol{\theta} \mid \boldsymbol{x}}}[-2\log f(\boldsymbol{x} \mid \boldsymbol{\theta})] + 2\left(\operatorname{E_{\boldsymbol{\theta} \mid \boldsymbol{x}}}[-2\operatorname{log}f(\boldsymbol{x} \mid \boldsymbol{\theta})] + 2\operatorname{log}f(\boldsymbol{x} \mid \operatorname{E_{\boldsymbol{\theta} \mid \boldsymbol{x}}}[\boldsymbol{\theta}])\right)$





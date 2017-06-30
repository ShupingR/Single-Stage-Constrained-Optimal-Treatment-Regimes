## Simulation Design for Normal and Mixture Normals

### Problem

maximize $Y=X^T\beta_1 + A X^T \beta_2 + \epsilon_Y$

minimize $Z=X^T\gamma_1 + A X^T \gamma_2 + \epsilon_Z$

Y is max when $A = sgn(X^T\beta_2)$

Z is min when $A = -sgn(X^T\gamma_2)$


### Factors under consideration and loss functions


#### 1. Probability of treatments disagree

Levels: 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9

$P = E[I\{(X^T\beta_2)(X^T\gamma_2)\}]$

$Loss_1 = (P - P.target)^2$


#### 2. Relative ratio of treatments regret

Levels: 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.75

$T = \frac{E{|X^T\beta_2|I\{(X^T\beta_2)(X^T\gamma_2)\}}/E|X^T\beta_2|}{E{|X^T\gamma_2| I\{(X^T\beta_2)(X^T\gamma_2)}/E{|X^T\gamma_2|}}$

$Loss_2 = (T - T.target)^2$


#### Solve to find the parameters $\beta_1$, $\beta_2$, $\gamma_1$ and $\gamma_2$

$Loss_{total} = Loss_1 + Loss_2$

Calculate total loss and solve by optim in R by methods SANN


#### 3. R-square

Level: 0.3, 0.5, 0.7


### Generate Data with covariates X from multivariate normal distribution

Generate $X \sim MVN(0, \Sigma)$, where
$$\mathbf{\Sigma} = \left[\begin{array}
{rr}
1 & 0 \\
0 & 1
\end{array}\right]$$
and add intercept 1 to $X$.


Generate $A$ from uniform random $\{-1, 1\}$.


Sample size $N = 10,000$.


#### Solve R-square analytically


$Y = X^T\beta_1 + A X^T \beta_2 + \epsilon_Y$
$= \beta_{10} + \beta_{11} X_1 + \beta_{12} X_2 + A\beta_{20} + A\beta_{21} X_1 + A\beta_{22} X_2 + \epsilon_Y$


$X \sim MVN(0, I)$


$Var(X_1) = Var(X_2) = 1$


$Var(A) = E(A^2) - E(A)^2 = 1$


$Var(AX_1) = E(A^2X_1^2) - [E(AX_1)]^2 = E(A^2X_1^2) - [E(A)E(X_1)]^2 = E(X_1^2) = Var(X_1) = 1$


Similarily, $Var(AX_2) = 1$


$SSR = Var(X^T\beta_1 + A X^T \beta_2)$
$= \beta_{11}^2Var(X_1) +  \beta_{12}^2Var(X_2) + \beta_{20}^2Var(A) + \beta_{21}^2Var(AX_1) + \beta_{22}^2Var(AX_2)$
$= \beta_{11}^2 + \beta_{12}^2 + \beta_{20}^2 + \beta_{21}^2 + \beta_{22}^2$


$SSE = Var(\epsilon) = \sigma_{Y}^2$


$SST = SSR + SSE$


$R^2 = SSR/SST = \frac{\beta_{11}^2 + \beta_{12}^2 + \beta_{20}^2 + \beta_{21}^2 + \beta_{22}^2}{\beta_{11}^2 + \beta_{12}^2 + \beta_{20}^2 + \beta_{21}^2 + \beta_{22}^2+ \sigma_{Y}^2} = SSR/(SSR + \sigma_Y^2)$


$\sigma_{Y}^2 = \frac{1 - R^2}{R^2} SSR = \frac{1 - R^2}{R^2} (\beta_{11}^2 + \beta_{12}^2 + \beta_{20}^2 + \beta_{21}^2 + \beta_{22}^2)$


### Generate Data with covariates X from multivariate mixture normal distributions


#### Solve R-square analytically


Recall $X \sim \sum_{j = 1}^K p_j N(\mu_j, \Sigma_j)$, where $X$ without intercept 1.


We write $X=Y+\epsilon$, where $Y = \mu_j$ with probability $p_j$, for $j=1,…,K$, and the conditional probability distribution of the $\epsilon$ given $Y$ is $N(0,\Sigma_j)$. Then we have


$E(X) = E(E(X\mid Y)) = E\left.\begin{cases} \vdots \\  \mu_j & \text{with probability } p_j \\  \vdots \end{cases}\right\} = \sum_{j=1}^K p_j\mu_j$


$Var(X) = E(Var(X \mid Y)) + Var(E(X \mid Y))$


$= E\left.\begin{cases} \vdots \\  \Sigma_j & \text{with probability }p_j \\  \vdots \end{cases}\right\} + Var\left.\begin{cases} \vdots \\ \mu_j & \text{with probability }p_j \\   \vdots \end{cases} \right\}$


$= \sum_{j=1}^K p_j \Sigma_j + \sum_{j=1}^K p_j(\mu_j-\bar\mu)(\mu_j-\bar\mu)^T$
, where $\bar\mu=\sum_{j=1}^K p_j\mu_j$




In our case,
$Y = X^T\beta_1 + A X^T \beta_2 + \epsilon_Y$
$= \beta_{10} + \beta_{11} X_1 + \beta_{12} X_2 + A\beta_{20} + A\beta_{21} X_1 + A\beta_{22} X_2 + \epsilon_Y$


Note $\bar{X} = (X_1, X_2)$ no intercept here and $X = (1, X_1, X_2)$ with intercept.


$\bar{X} \sim \sum_{j=1}^{K} p_{j}BVM(\mu_j, \Sigma_j)$ and $\bar{\beta_1} = (\beta_{11}, \beta_{12})$


Then, $\bar{X}$ has mean $\mu_{\bar{X}} = E{\bar{X}} = \sum_{j=1}^K p_j\mu_j$ and covariance matrix $\Sigma_{\bar{X}} = \sum_{j=1}^K p_j \Sigma_j + \sum_{j=1}^K p_j(\mu_j-\bar\mu)(\mu_j-\bar\mu)^T$, where $\bar\mu=\sum_{j=1}^K p_j\mu_j$


ie $\bar{X} \sim Dist'n(\mu_{\bar{X}}, \Sigma_{\bar{X}})$


For our potential outcome random scalar variable$Y$ and $X$ representing a random vector from mixture gaussian model, then we have


$Y = X^T\beta_1 + A X^T \beta_2 + \epsilon_Y$
$= \beta_{10} + \beta_{11} X_1 + \beta_{12} X_2 + A\beta_{20} + A\beta_{21} X_1 + A\beta_{22} X_2 + \epsilon_Y$

and let $\bar{Y} =  X^T\beta_1 + A X^T \beta_2$, then $SSR = Var(\bar{Y}) = E(\bar{Y}^2) -E(\bar{Y})^2$


As $A$ independent of $X$, $E(A)=0$, $A^2=1$,

and also $Var(\bar{X}) = E(\bar{X}\bar{X}^T) - E(\bar{X})E(\bar{X})^T$, $E(\bar{X}\bar{X}^T) = Var(\bar{X}) + E(\bar{X})E(\bar{X})^T = \Sigma_{\bar{X}} + \mu_{\bar{X}}\mu_{\bar{X}}^T$


$\bar{Y}^2 = (\beta_1^{T}X + A\beta_2^{T}X)^2 = (\beta_1^{T}X)^2 + 2A(\beta_1^{T}X)(\beta_2^{T}X) + A^2(\beta_1^{T}X)^2$


$E(\bar{Y}^2) = E\{(\beta_1^{T}X)^2 + 2A(\beta_1^{T}X)(\beta_2^{T}X) + A^2(\beta_1^{T}X)^2\} = E\{(\beta_1^{T}X)^2\} + 0 + E\{(\beta_2^{T}X)^2\}$


$= E\{ \beta_1^T X X^T \beta_1\} +E\{ \beta_2^T X X^T \beta_2\} = E \{(\beta_1 +\beta_2)^T X X^T (\beta_1 + \beta_2) \} = (\beta_1 +\beta_2)^T E (X X^T) (\beta_1 + \beta_2)$,


where  
$${E (X X^T)} = \left[\begin{array}
{rrr}
1 & \mu_{X_1} & \mu_{X_2} \\
\mu_{X_1} &  &  \\
\mu_{X_2} &  & (\Sigma_{\bar{X}} + \mu_{\bar{X}}\mu_{\bar{X}}^T)
\end{array}\right]
$$


$[E(\bar{Y})]^2 = [E(\beta_1^{T}X + A\beta_2^{T}X)]^2 = [\beta_1^T E(X)]^2 = (\beta_{10} + \bar{\beta_1}\mu_{\bar{X}})^2$


Therefore,

$SSR = Var(\bar{Y}) = E(\bar{Y}^2) - [E(\bar{Y})]^2 =  (\beta_1 +\beta_2)^T E (X X^T) (\beta_1 + \beta_2) - (\beta_{10} + \bar{\beta_1}\mu_{\bar{X}})^2$


$SST = Var(\bar{Y}) + \sigma^2_{\epsilon_Y}$


$R^2 = SSR/SST = SSR/ (SSR + \sigma^2_{\epsilon_Y})$


$\sigma_{Y}^2 = \frac{1 - R^2}{R^2} SSR = \frac{1 - R^2}{R^2} \{(\beta_1 +\beta_2)^T E (X X^T) (\beta_1 + \beta_2) - (\beta_{10} + \bar{\beta_1}\mu_{\bar{X}})^2 \}$


### Solving for Multivariate normal X by SANN

```{r results="asis",echo=FALSE, id:"iz65h16i"}
rm(list=ls())
library(xtable)
dir <- "/home/sruan/Workspace/Research/constrained-optimal-regime-1/simulation-design/"
setwd(dir)
parm <- read.table("parmSolList_SANN.txt")
tab1 <- xtable(parm, caption="Parameters Solution by SANN",  digits=rep(4,15),
              caption.placement = "top")
print(tab1, type="html")
```

<br>

```{r results="asis",echo=FALSE, id:"iz65h16x"}
loss <- read.table("lossList_SANN")
tab2 <- xtable(loss, caption="Total Loss by SANN",  digits=c(0, 5, 5, 5),
              caption.placement = "top")
print(tab2, type="html")
```

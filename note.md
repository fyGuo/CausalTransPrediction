\documentclass[11pt]{article}

\usepackage[margin=1in]{geometry}
\usepackage{amsmath, amssymb, amsthm}
\usepackage{tikz}
\usetikzlibrary{arrows.meta, bending}

\usepackage{setspace}
\usepackage{graphicx}
\usepackage{caption}

\newcommand{\E}{\mathbb{E}}
\newcommand{\indep}{\perp\!\!\!\perp}
\newcommand{\pp}{\mathbb{P}}

\newtheorem{theorem}{Theorem}
\newtheorem{corollary}{Corollary}
\title{Notes}
\author{Author Name}
\date{\today}

\begin{document}

\maketitle

\section*{Counterfactual, transportable risk prediction in a time-fixed setting}

Let us work with a simple case with time-fixed data. Following the notation used in the main text, the time-fixed data is iid and defined $O_i\equiv (X_i, P_i, L_i, S_i,A_i, Y_i)$ where $X_i$ is a set of baseline covariates, $P_i$ is the risk predictors and $P_i \subset X_i$, $L_i$ is confounders, $S_i$ is the cohort indicator, with $\mathcal{S}_i = 0, 1, \ldots, K$. We use $S_i = 0$ to indicate that people are in the target population. Exposures and outcomes are denoted by $A_i$ and $Y_i$. For notational simplicity, we omit the subscript and use the generic notation. For $k = 0,\ldots, K$, we select a subset $X^{(k)} \subset X$ measured in both $k$-th cohort and target population to ensure transportability. We also select a subset $L^{(k)} \in L$ measured in the $k$-th cohort to ensure conditional randomization. The causal estimand is a function $\E(Y^{a}\mid P, S=0)$. Below, we show how to identify by combining the $k$-th cohort with the target population. We need the following identification assumptions.


\begin{enumerate}

\item[\textbf{A1.}] \label{a:consistency}\textbf{Consistency.}
If $A_i = a$, then $Y^a_i = Y_i$. Database participation does not itself affect the outcome.

\item[\textbf{A2.}] \label{a:exchange0}\textbf{No unmeasured baseline confounding in cohort $k$}
For each database $k$,
\begin{equation}
Y^a \;\indep\; A \;\mid\; L^{(k)}, X^{(k)}, S=k
\label{eq:exchange0}
\end{equation}
Within database $k$, ${L}^{(k)}$ and ${X}^{(k)}$ are sufficient to render baseline treatment conditionally independent of the counterfactual outcome.

\textcolor{red}{Fuyu: this assumption is different from the original one, we can modify it if it is two strong. I think the most important thing is that $X$ cannot be a mediator or collider from $A$ to $Y$}

\item[\textbf{A3.}] \label{a:pos_treat0}\textbf{Positivity of treatment (time-fixed).}
For each database $k$,
\begin{equation*}
\Pr\!\bigl[A = a \;\big|\; L^{(k)} = l,\; X^{(k)} = x, S = k\bigr] > 0
\end{equation*}
for all $a, l, x$ in the support.

\item[\textbf{A4.}] \label{a:transport}\textbf{Transportability.}
For each $k$, combing the target population and the cohort $k$, it follows
\begin{equation}
Y^a \;\indep\; S \;\mid \; X^{(k)}.
\label{eq:transport}
\end{equation}
Alternatively, it can be expressed as $\Pr(S=k\mid X^{(k)}, Y^a) = \Pr(S=0 \mid X^{(k)}, Y^a)$


\item[\textbf{A5.}] \label{a:pos_study}\textbf{Positivity of database membership.}
For each $k$, 
\begin{equation*}
0 < \Pr\!\bigl[S = k \;\big|\; \bmath{X}^{(k)} = x^{(k)}\bigr] < 1
\end{equation*}
for every $x^{(k)}$ in the support of the target population.

\end{enumerate}
One DAG satisfying the above assumptions are shown below
\begin{figure}[htbp]
\centering
\begin{tikzpicture}[> = stealth, shorten > = 1pt, auto, node distance = 1.75cm, semithick]
\tikzstyle{every state}=[
  draw = black,
  thick,
  fill = white
]
\node (x) {$X$};
\node (s) [right of=x] {$S$};
\node (a) [right of=s] {$A$};
\node (y) [right of=a] {$Y$};
\node (l) [below of=x] {$L$};

% Main horizontal path
\path[->] (x) edge node {} (s);
\path[->] (s) edge node {} (a);
\path[->] (a) edge node {} (y);

% X -> L, L -> A, L -> Y
\path[->] (x) edge node {} (l);
\path[->] (l) edge node {} (a);
\path[->] (l) edge node {} (y);

% X curved top arcs to A and Y
\path[->] (x) edge [out=40, in=-210] node {} (a);
\path[->] (x) edge [out=40, in=-210] node {} (y);

\end{tikzpicture}
\caption{A causal directed acyclic graph relating $X$, $S$, $L$, $A$, and $Y$.}
\label{fig:dag-timefixed}
\end{figure}

\caption{A causal directed acyclic graph relating $X$, $S$, $A$, $Y$, and $L$.}
\label{fig:dag}

Now, we are able to give the identification results.

\begin{theorem}
    Under Assumptions A1 to A5, it follows that
    \begin{align*}
        \E(Y^a\mid P, S=0) = \E\left[\E\left\{\E\left(Y\mid A=a, X^{(k)}, L^{(k)}, S=k\right)\mid X^{(k)}, S=0\right\}\mid P, S=0\right]
    \end{align*}
\end{theorem}
\begin{proof}
    \begin{align*}
        &\E(Y^a\mid P, S=0)\\
        &=\E\left\{\E(Y^a\mid P, X^{(k)}, S=0)\mid P, S=0\right\}\\
        &=\E\left\{\E(Y^a\mid  X^{(k)}, S=0)\mid P, S=0\right\}&& P\subset X^{(k)} \\
        &= \E\left\{\E(Y^a\mid  X^{(k)}, S=k)\mid P, S=0\right\}&& \textrm{Assumption A5}\\
        &= \E\left[\E\left\{\E(Y^a\mid  X^{(k)}, S=k, L^{(k)})\mid  X^{(k)}, S=k\right\}\mid P, S=0\right]\\
        &= \E\left[\E\left\{\E(Y^a\mid  X^{(k)}, S=k, L^{(k)}, A=a)\mid  X^{(k)}, S=k\right\}\mid P, S=0\right] && \textrm{Assumption A2}\\
        &= \E\left[\E\left\{\E(Y\mid  X^{(k)}, S=k, L^{(k)}, A=a)\mid  X^{(k)}, S=k\right\}\mid P, S=0\right] && \textrm{Assumption A1}
    \end{align*}
The positivity assumptions ensure the conditional expectations are well-defined.
\end{proof}
The above Theorem implies an estimation algorithm.
\begin{enumerate}
    \item In study $k$, build a model to estimate $\E(Y\mid  X^{(k)}, S=k, L^{(k)}, A=a)$, denoted by $\mu^{(k)}(X^{(k)}, L^{(k)})$ (omit small$a$ as it is a constant).
    \item For each individual $j$ in cohort $k$ fit that model and get $\hat \mu^{(k)}_j(X^{(k)}_j, L^{(k)}_j)$.
    \item Fit a model, say regression, for in cohort $k$ by regressing $\mu^{(k)}_j(X^{(k)}_j, L^{(k)}_j)$ on $X^{(k)}_j$, and get a model denoted by $\tilde \mu^{(k)}(X^{(k)})$.
    \item In the target population, for each subject $i$, plug in $X^{(k)}_i$ and get $\tilde \mu^{(k)}_i (X^{(k)}_i)$.
    \item Fit a model by regressing $\tilde \mu^{(k)}_i (X^{(k)}_i)$ on $P_i$ in the target population, and get all the coefficients.
\end{enumerate}
If people want, they can meta-analyze the coefficient from $K$ cohorts.

\textcolor{red}{Sinclair: let's deal with the meta-analysis case later}

\newpage
\section{Linear models and asymptotic variance}
To do the analysis, we usually make a linear model for $\E(Y^0\mid P, S=0)$, 
$$\E(Y^0\mid P, S=0) = \beta^\top P ,$$
where $\beta \in \mathbb{R}^p$. We would like to estimate and make an inference on $\beta^\top$. First, we combine each cohort with the target population and estimate $\hat \beta$. We then use meta-analysis to combine $\hat \beta$ across $K$ estimates. In the first stage, we combine cohort $k$ with the target population, so the data is defined as $O\equiv (P, X^{(k)}, L^{(k)}, S^{(k)}, A, Y)$, where $S^{(k)} = 1$ if the data comes from the cohort $k$ and 0 for the target population. The following theorem gives the set of influence functions of $\beta$.

\begin{theorem}\label{theorem:if_time_fixed}
Under the condition that $\E\left\{(1-S^{(k)} P P^\top)\right\}$ is invertible. The set of influence functions for $\beta$ is 

$$\{\phi_\beta(O;d): d(P) \text{  is an arbitrary function of }P\},$$

where 
\begin{align*}
 \phi_\beta(O)&=V^*(d)^{-1}
d(P)\Bigg[
(1-S^{(k)})\left\{\psi^*(X^{(k)})-\beta^{*\top}P\right\} \\
& \quad +
\frac{S^{(k)}\rho^*_0(X^{(k)})}{\rho^*_k(X^{(k)})}
\left\{
\frac{1-A}{1-\pi^*(X^{(k)},L^{(k)})}\{Y-m^*(X^{(k)},L^{(k)})\} +m^*(X^{(k)},L^{(k)})- \psi^*(X^{(k)})
\right\}
\Bigg],\\
V(d) &= \E\left\{(1-S^{(k)} d(P) P^\top)\right\},\\
m(X^{(k)}, L^{(k)})&= \E(Y \mid A, X^{(k)}, L^{(k)}, S_k = 1),\\
\psi(X^{(k)}) &= \E\left\{m(X^{(k)}, L^{(k)}) \mid X, S_k = 1\right\},\\
\rho_0(X^{(k)}) &= \E(S^{(k)} = 0\mid X^{(k)}),\\
\rho_k(X^{(k)}) &= \E(S^{(k)} = 1\mid X^{(k)}),\\
\pi(X^{(k)}, L^{(k)}) &= \Pr(A=1\mid X^{(k)}, L^{(k)}, S^{(k)} = 1)
\end{align*}
\end{theorem}


Before giving the proof, let us look at how different people contribute to the estimation.
\begin{enumerate}
    \item Among people in the target population , $S^{(k)} = 0$, the contribution to the unscaled influence function is $\psi^*(X^{(k)})-\beta^{*\top}P$.
    \item Among people in the cohort $k$,  the contribution is $$\frac{S^{(k)}\rho^*_0(X^{(k)})}{\rho^*_k(X^{(k)})}
\left\{
\frac{1-A}{1-\pi^*(X^{(k)},L^{(k)})}\{Y-m^*(X^{(k)},L^{(k)})\} +m^*(X^{(k)},L^{(k)})- \psi^*(X^{(k)})\right\},$$
a reweighted version of AIPW estimator.The weight is $\frac{S^{(k)}\rho^*_0(X^{(k)})}{\rho^*_k(X^{(k)})}$. It is very similar to the policy weighting in the reinforcement literature but replacing treatment with the cohort indicator.
\end{enumerate}
It can be seen that cohort $k$ people contribute to the AIPW estimation because they have outcome measured, while target population contribute to the $\beta$ estimation. The function $\psi^*(X^{(k)})$ bridges the two populations. It follows that

\begin{align*}
    \psi^*(X^{(k)}) &= \E\left\{m^*(X^{(k)}, L^{(k)})\mid X^{(k)}, S^{k} = 1\right\}\\
    &=  \E\left[\E\left\{Y\mid S^{(k)}= 1, A= 0, X^{(k)}, L^{(k)}\right\}\mid X^{(k)}, S^{k} = 1\right]\\
    &=  \E\left[\E\left\{Y^0\mid S^{(k)}= 1, X^{(k)}, L^{(k)}\right\}\mid X^{(k)}, S^{k} = 1\right]\\
    &= \E\left(Y^0\mid X^{(k)}, S^{k} = 1\right)\\
    & = \E\left(Y^0\mid X^{(k)}\right),
\end{align*}
which indicates that although $\psi^*(X^{(k)})$ is defined among people with $S^{(k)} = 1$, it does not actually depend on $S^{(k)}$. This property is the foundation of the transportation between populations.
\begin{proof}
    According to the linear model, it follows that
    \begin{align*}
        \E\left[d(P) \left\{\beta^{*\top} P- \mu^*_0(P)\right\}\right] = 0
    \end{align*}
    for arbitrary functions $d$ of $P$. It follows that
    \begin{align*}
        \E\left(d(P) \left[\beta^{*\top} P- \E\left\{\psi(X^{(k)})\mid P, S=0\right\}\right]\right) = 0,
    \end{align*}
    and 
     \begin{align*}
        \E\left(d(P)(1-S^{(k)}) \left[\beta^{*\top} P- \E\left\{\psi(X^{(k)})\mid P, S^{(k)}\right\}\right]\right) = 0,
    \end{align*}
    which leads to
      \begin{align*}
        \E\left(d(P)(1-S^{(k)}) \left[\beta^{*\top} P- \psi^*(X^{(k)})\right]\right) = 0,
    \end{align*}
    Solving the equation it follows that,
    $\beta^{*} = V(d)^{-1} \E\left\{(1-S^{(k)})d(P)\psi(X^{(k)})\right\}$. Let us define $\theta \equiv \E\left\{(1-S^{(k)})d(P)\psi^*(X^{(k)})\right\}$.
    Using the chain rule of influence functions,
    $\phi_\beta(O) = V(d)^{-1} \phi_\theta(O) + V(d)^{-1}\phi_V(O)V(d)^{-1}\theta$. Let us derive the influence function for $\theta$ first.
    \begin{align*}
        \theta &= \E\left\{(1-S^{(k)})d(P)\psi^*(X^{(k)})\right\}\\
        &=  \E\left\{\rho_0(X^{(k)})d(P)\psi^*(X^{(k)})\right\}\\
        &= \E\left[d(P)\rho_0(X^{(k)})\E\left\{m(X^{(k)}, L^{(k)})\mid X^{(k)}, S_k = 1\right\} \right]\\
        &= \E\left[d(P)\rho_0(X^{(k)})\E\left\{\frac{S^{(k)}}{\rho_k(X^{(k)})}m(X^{(k)}, L^{(k)})\mid X^{(k)}\right\} \right]\\
        &= \E\left\{d(P)\rho_0(X^{(k)})\frac{S^{(k)}}{\rho_k(X^{(k)})}m(X^{(k)}, L^{(k)})\right\}\\
        & = \E\left\{d(P)\rho_0(X^{(k)})\frac{S^{(k)}}{\rho_k(X^{(k)})}\E(Y\mid A=0, X^{(k)}, L^{(k)}, S^{(k)})\right\}\\
         & = \E\left[d(P)\rho_0(X^{(k)})\frac{S^{(k)}}{\rho_k(X^{(k)})}\E\left\{\frac{1-A}{1-\pi(X^{(k)}, L^{(k)})}Y\mid  X^{(k)}, L^{(k)}, S^{(k)}\right\}\right]\\
         & = \E\left\{d(P)\rho_0(X^{(k)})\frac{S^{(k)}}{\rho_k(X^{(k)})}\frac{1-A}{1-\pi(X^{(k)}, L^{(k)})}Y\right\}\\
         & = \E\left\{d(P)S^{(k)}\frac{1-\rho_k(X^{(k)})}{\rho_k(X^{(k)})}\frac{1-A}{1-\pi(X^{(k)}, L^{(k)})}Y\right\}
    \end{align*}
    Using the perturbation technique, 


    \begin{align*}
        \frac{\partial }{\partial t}\theta_t|_{t=0} &=\left.\frac{\partial}{\partial t} \E_t\left\{d(P)S^{(k)}\frac{1-\rho^*_k(X^{(k)})}{\rho^*_k(X^{(k)})}\frac{1-A}{1-\pi^*(X^{(k)}, L^{(k)})}Y\right\}\right|_{t=0}\\
        &\quad + \E\left\{d(P)S^{(k)}\left.\frac{\partial }{\partial t}\frac{1-\rho_k(X^{(k)})}{\rho_k(X^{(k)})}\right|_{t=0}\frac{1-A}{1-\pi^*(X^{(k)}, L^{(k)})}Y\right\}\\
        &\quad + \E\left\{d(P)S^{(k)}\frac{1-\rho^*_k(X^{(k)})}{\rho^*_k(X^{(k)})}\left.\frac{\partial }{\partial t}\frac{1-A}{1-\pi(X^{(k)}, L^{(k)})}\right|_{t=0}Y\right\}\\
        &= (i) + (ii) + (iii)
    \end{align*}

Let is use $\mathcal{S}(\cdot)$ to denote the score function. For term (i),

$$(i) = \E\left[\left\{d(P)S^{(k)}\frac{1-\rho^*_k(X^{(k)})}{\rho^*_k(X^{(k)})}\frac{1-A}{\pi^*(X^{(k)}, L^{(k)})}Y-\theta\right\}\mathcal{S}(O)\right]$$
For term (ii) 
    \begin{align*}
        (ii) &= -\E\left\{d(P)S^{(k)}\frac{1}{\rho^{*2}_k(X^{(k)})}\frac{\partial }{\partial t}\E_t(S^{(k)}\mid X^{(k)})|_{t=0}\frac{1-A}{\pi^*(X^{(k)}, L^{(k)})}Y\right\}\\
     &= -\E\left[d(P)S^{(k)}\frac{1}{\rho^{*2}_k(X^{(k)})}\frac{\partial }{\partial t}\E\left\{S^{(k)}\mathcal{S}(S^{(k)}\mid X^{(k)})\mid X^{(k)}\right\}|_{t=0}\frac{1-A}{\pi^*(X^{(k)}, L^{(k)})}Y\right]\\
     & = -\E\left[d(P)\frac{1}{\rho^{*}_k(X^{(k)})}S^{(k)}\mathcal{S}(S^{(k)}\mid X^{(k)})\psi^*(X^{(k)})\right]\\
     &= -\E\left[d(P)\frac{1}{\rho^{*}_k(X^{(k)})}S^{(k)}\E\left\{\mathcal{S}(O)\mid X^{(k)},{S^{(k)}}\right\}\psi^*(X^{(k)})\right]\\
     &\quad +  \E\left[d(P)\frac{1}{\rho^{*}_k(X^{(k)})}S^{(k)}\E\left\{\mathcal{S}(O)\mid X^{(k)}\right\}\psi^*(X^{(k)})\right]\\
     &=  -\E\left[d(P)\frac{1}{\rho^{*}_k(X^{(k)})}\left\{S^{(k)} -\rho^{*}_k(X^{(k)}) \right\}
     \psi^*(X^{(k)})\mathcal{S}(O)\right].
    \end{align*}
For term (iii), one can use the same trick for term (ii) and show that 
$$(iii) = \E\left[d(P) \frac{1-\rho_k^*(X^{(k)})}{\rho_k^*(X^{(k)})}S^{(k)}\frac{A-\pi(X^{(k)}, L^{(k)})}{1- \pi(X^{(k)}, L^{(k)})}m (X^{(k)}, L^{(k)})\mathcal{S}(O)\right].$$
Summing up (i), (ii), and (iii), with some algebra, it follows that

\begin{align*}
    \phi_\theta(O)
&=d(P)\left(
(1-S^{(k)})\psi(X^{(k)}) - \theta \right.\\
&\left. \quad +\frac{S^{(k)}\rho_0(X^{(k)})}{\rho_k(X^{(k)})}
\left[
\frac{1-A}{1-\pi(X^{(k)}, L^{(k)})}\{Y-m(X^{(k)}, L^{(k)})\}+m(X^{(k)}, L^{(k)}) - \psi(X^{(k)}, L^{(k)})
\right]\right)
\end{align*}

Now, let us look at the influence function for $V(d)$, 
\begin{align*}
    \frac{\partial }{\partial t} V(d)|_{t=0} &= \frac{\partial }{\partial t}\E_t\left\{(1-S^{(k)})d(P)P^\top\right\}|_{t=0}\\
    &=\E\left\{(1-S^{(k)})d(P)P^\top\mathcal{S}(O)\right\}
\end{align*}
Thus, 
$\phi_V(O) =(1-S^{(k)})d(P)P^\top - V^*(d) $
Now using the chain rule, it follows that
\begin{align*}
    \phi_\beta(O) &=V^*(d)^{-1} \phi_\theta(O) + V^*(d)^{-1}\phi_V(O)V^*(d)^{-1}\theta^*\\
    &= V^*(d)^{-1} \phi_\theta(O) + V^*(d)^{-1}\phi_V(O)\beta^*\\
    &= V^*(d)^{-1}\left\{\phi_\theta(O) + \phi_V(O)\beta^*\right\}\\
    &= V^*(d)^{-1}d(P)\left\{\phi_\theta(O) + \phi_V(O)\beta^*\right\}\\
    &= V^*(d)^{-1}
d(P)\BiggV^*(d)^{-1}
d(P)\Bigg[
(1-S^{(k)})\left\{\psi^*(X^{(k)})-\beta^{*\top}P\right\} \\
& \quad +
\frac{S^{(k)}\rho^*_0(X^{(k)})}{\rho^*_k(X^{(k)})}
\left\{
\frac{1-A}{1-\pi^*(X^{(k)},L^{(k)})}\{Y-m^*(X^{(k)},L^{(k)})\} +m^*(X^{(k)},L^{(k)})- \psi^*(X^{(k)})
\right\}
\Bigg]
\end{align*}
\end{proof}

Based on the influence functions, we are ready to give its double robustness property.
\begin{corollary}
    Under the conditions that
    \begin{enumerate}
        \item Either $m^*(X^{(k)}, L^{(k)})$ or $\pi^*(X^{(k)}, L^{(k)})$ is correctly specified.
        \item and Either $\psi^*(X^{(k)}, L^{(k)})$ or $\rho_0^*(X^{(k)}, L^{(k)})$ is correctly specified.
    \end{enumerate}
    The estimator $\hat \beta$ that solves the estimating equation 
    $\pp_n \phi_\beta(O_i) = 0$ is unbiased.
\end{corollary}
\begin{proof}
    It suffices to show that 
    $\E\phi_\beta(O)$ is unbiased under the conditions. Under condition 1, it is well known that
    \begin{align*}
        \E\left[\frac{1-A}{1-\pi(X^{(k)},L^{(k)})}\{Y-m(X^{(k)},L^{(k)})\} + m(X^{(k)},L^{(k)})\mid X^{(k)}, L^{(k)}\right] = m^*(X^{(k)},L^{(k)})
    \end{align*}
    Now suppose $\psi(X^{(k)}) = \psi^*(X^{(k)})$, it follows that 
    \begin{align*}
        \E\left\{\phi_\beta(O)\right\} &= \left(V^*(d)^{-1}
d(P)\E\left[(1-S^{(k)})\left\{\psi^*(X^{(k)})-\beta^{*\top}P\right\}\right]\right)\\
&\quad  +\E\left(V^*(d)^{-1}
d(P) \frac{S^{(k)}\rho_0(X^{(k)})}{\rho_k(X^{(k)})}
\left[\E\{m^*(X^{(k)},L^{(k)})\mid X^{(k)}, S^{(k)} =1\}- \psi^*(X^{(k)})
\right]\right)\\
&= 0
    \end{align*}
Now suppose $\rho_0(X^{(k)}) = \rho^*_0(X^{(k)}) $, it follow that 
\begin{align*}
        \E\left\{\phi_\beta(O)\right\} &=\E\left( V^*(d)^{-1}
d(P)\E\left[(1-S^{(k)})\left\{\psi(X^{(k)}) - \beta^{*\top}P\right\}\right]\right)\\
&\quad  +\E\left(V^*(d)^{-1}
d(P) \frac{S^{(k)}\rho^*_0(X^{(k)})}{\rho^*_k(X^{(k)})}
\left[\E\{m^*(X^{(k)},L^{(k)})\mid X^{(k)}, S^{(k)} =1\}- \psi(X^{(k)})
\right]\right)\\
&=  \E\left(V^*(d)^{-1}
d(P)\E\left[(1-S^{(k)})\left\{\psi(X^{(k)})-\beta^{*\top}P\right\}\right]\right)\\
&\quad  +\E\left[V^*(d)^{-1}
d(P) \rho^*_0(X^{(k)})
\left\{\psi^*(X^{(k)})- \psi(X^{(k)})
\right\}\right]\\
&=  \E\left(V^*(d)^{-1}
d(P)\E\left[\rho^*_0(X^{(k)})\left\{\psi(X^{(k)})-\beta^{*\top}P + \psi^*(X^{(k)})- \psi(X^{(k)})\right\}\right]\right)\\
&=\E\left(V^*(d)^{-1}
d(P)\E\left[\rho^*_0(X^{(k)})\left\{-\beta^{*\top}P + \psi^*(X^{(k)})\right\}\right]\right)\\
&=\E\left(V^*(d)^{-1}
d(P)\E\left[(1-S^{(k)})\left\{-\beta^{*\top}P + \psi^*(X^{(k)})\right\}\right]\right)\\
&=\E\left(V^*(d)^{-1}
d(P)(1-S^{(k)})\E\left[\left\{-\beta^{*\top}P + \psi^*(X^{(k)})\right\}\mid  S^{(k)} = 0, P\right]\right)\\
&= 0
    \end{align*}

\end{proof}


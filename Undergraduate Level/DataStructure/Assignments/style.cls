%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% style class for LaTeX files               v0.9 March 18 2024   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{style}[2025/01/01 v0.9 style class for LaTeX files]
\LoadClass[12pt,a4paper]{article}

% --- PACKAGES FROM preamble.tex ---
% Figures, math, tables, fonts, references, etc.
\RequirePackage{subfigure,graphicx}
\RequirePackage{amsmath,amsfonts,amsthm,latexsym,amssymb,euscript,xr}
\RequirePackage{booktabs}
\RequirePackage[nodayofweek]{datetime}
\RequirePackage{hyperref}
\RequirePackage[version=4]{mhchem}
%\RequirePackage[protrusion=true,expansion=true]{microtype}
\RequirePackage{siunitx}
    \sisetup{detect-weight,detect-mode}
\RequirePackage{wrapfig}
    \setlength\intextsep{2pt}
\RequirePackage[table]{xcolor}
\RequirePackage{color,colortbl,tabularx}
\RequirePackage{pifont}
\RequirePackage[numbib]{tocbibind}
    \renewcommand{\tocbibname}{References relevant to the experiment description}
\RequirePackage{float}
\RequirePackage[utf8]{inputenc}
\RequirePackage{pgf,tikz}
\RequirePackage{tikz-3dplot}
\RequirePackage{multirow}
\RequirePackage{setspace}
\RequirePackage{bigstrut}
\RequirePackage{etoolbox}
\RequirePackage{cancel}
\RequirePackage{ragged2e}
\RequirePackage{xparse}
\RequirePackage{bm}
\RequirePackage{boxedminipage}
\RequirePackage{listings}
\RequirePackage{ifthen}

% TikZ libraries and layers
\usetikzlibrary{calc,matrix,patterns,shadings,backgrounds,intersections,decorations.pathreplacing,math}
\pgfdeclarelayer{myback}
\pgfdeclarelayer{background}
\pgfsetlayers{myback,background,main}
\tikzset{myfillcolor/.style = {draw,fill=#1}}
\NewDocumentCommand{\fillpattern}{O{north west lines} O{blue!50} m m}{%
\draw[pattern=#1, pattern color=#2] (#3.north west)rectangle (#4.south east);
}

% --- FONT ---
\usepackage{helvet}
\renewcommand{\familydefault}{\sfdefault}

% --- COLOR DEFINITIONS ---
\definecolor{LightColor}{rgb}{0.8,0.77,0.9}
\def\CC{\cellcolor{red!30}}
\def\CB{\cellcolor{green!30}}
\def\RC{\color{red}}
\def\BC{\color{blue}}
\def\GC{\color{green}}
\renewcommand{\CancelColor}{\color{red}}

% --- SPECIAL COMMANDS ---
\def\code#1{\texttt{#1}}
\newcommand{\BOX}[2]{\fbox{\begin{minipage}{#1}\centering #2 \end{minipage}}}
\newcommand{\WIDEBOX}[1]{\vspace*{0.5cm}\fbox{\begin{minipage}{0.99\textwidth} #1 \end{minipage}}\vspace*{0.5cm}}
\newcommand\FLine[4]{( ( (#3) - (#1)*(#4) )/(#2) )}
\newcommand\InvFLine[4]{(( (#3)-(#2)*(#4))/((#1)))}
\newcommand{\DrawLine}[8]{%
    \pgfmathparse{#2  == 0.0}
        \ifnum\pgfmathresult=0 
            \coordinate (A) at ({#4}  , {\FLine{#1}{#2}{#3}{ #4} });
        \coordinate (B) at ({#5}  ,{\FLine{#1}{#2}{#3}{ #5 }});
    \else
        \coordinate (A) at (#3/#1,#4);
        \coordinate (B) at (#3/#1,#5);
        \fi
        \draw[#8] (A) -- (B) node[right]{ #7} coordinate[midway](M);
    \pgfmathparse{#6  == 0.0}
        \ifnum\pgfmathresult=0 
            \draw [#8,thick,solid,->] (M) -- ($(M)! #6cm!90:(A)$);
        \fi
}
\newcommand{\DisplayCoordinates}[6]{%
    \pgfmathparse{#1  == 0.0}
    \ifnum\pgfmathresult=1 
            \draw[#6] (#5,#2) -- (#1,#2) -- node[#4]{#2} (-#5,#2);
    \else
        \pgfmathparse{#2  == 0.0}
        \ifnum\pgfmathresult=1
            \draw[#6] (#1,#5) -- (#1,#2) -- node[#3]{#1} (#1,-#5);
        \else
            \draw[#6] (#1,#2) -- (#1,0) -- node[#3]{#1} (#1,-#5);
            \draw[#6] (#1,#2) -- (0,#2) -- node[#4]{#2} (-#5,#2);
        \fi
        \fi
}
\newcommand{\DisplayX}[4]{%
    \DisplayCoordinates{#1}{0}{#2}{}{#3}{#4};
}
\newcommand{\DisplayY}[4]{%
    \DisplayCoordinates{0}{#1}{}{#2}{#3}{#4};
}

% --- ENVIRONMENTS ---
\newenvironment{practice}
    {\par\medskip
     \noindent\begin{boxedminipage}{0.98\textwidth}\textbf{Practice Exercise:}\par\medskip}
    {\end{boxedminipage}\medskip}

\newenvironment{instructions}
    {\par\medskip
     \noindent\begin{boxedminipage}{0.98\textwidth}\textbf{Instructions:}\par\medskip}
    {\end{boxedminipage}\medskip}

\NewDocumentEnvironment{example}{o}{%
  \par\medskip
  \noindent\begin{boxedminipage}{0.98\textwidth}%
  \textbf{Example%
    \IfValueT{#1}{: #1}%
  }\par\medskip
}{%
  \end{boxedminipage}\medskip
}

\NewDocumentEnvironment{testcase}{o}{%
  \par\medskip
  \noindent\begin{boxedminipage}{0.98\textwidth}%
  \textbf{Test Case%
    \IfValueT{#1}{: #1}%
  }\par\medskip
}{%
  \end{boxedminipage}\medskip
}

% --- LISTINGS STYLE ---
\lstdefinestyle{customcpp}{
    language=C++,
    frame=single,
    basicstyle=\ttfamily\small,
    keywordstyle=\color{blue},
    commentstyle=\color{gray},
    stringstyle=\color{red},
    breaklines=true,
    breakatwhitespace=true
}
\lstnewenvironment{codelisting}
  { \lstset{style=customcpp} }
  {}

% --- SECTIONING, TITLE, LAYOUT (from original style.cls) ---
\renewcommand{\@maketitle}{
 \thispagestyle{empty}
 \begin{flushleft}%
  {
   \fontsize{12}{12}\selectfont
   \@title
   \par
  }
 \end{flushleft}
}

\def\author{\@ifnextchar[{\@getshortauthor}{\@getauthor}}
\def\@getauthor#1{\gdef\@author{#1}\gdef\@shortauthor{#1}}
\def\@getshortauthor[#1]#2{\gdef\@shortauthor{#1}\gdef\@author{#2}}
\newcommand{\address}[1]{\gdef\@address{#1}}
\renewcommand{\abstract}[1]{\gdef\@abstract{#1}}
\newcommand{\heading}[1]{\gdef\@heading{#1}}

\def\sectionLevel{1}
\def\sectionIndent{0mm}
\def\sectionSpaceBefore{12pt}
\def\sectionSpaceAfter{1pt}
\def\sectionStyle{\normalsize\bfseries}
\renewcommand{\section}{%
  \@startsection{section}{\sectionLevel}{\sectionIndent}{\sectionSpaceBefore}%
                {\sectionSpaceAfter}{\sectionStyle}}%
\renewcommand{\thesection}{\Alph{section}}

\def\sectionLevelsub{2}
\def\sectionIndentsub{0mm}
\def\sectionSpaceBeforesub{12pt}
\def\sectionSpaceAftersub{1pt}
\def\sectionStylesub{\normalsize}
\renewcommand{\subsection}{%
  \@startsection{subsection}{\sectionLevelsub}{\sectionIndentsub}{\sectionSpaceBeforesub}%
                {\sectionSpaceAftersub}{\sectionStylesub}}%

\renewcommand{\@makecaption}[2]{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{\normalsize\fontsize{10}{12}\selectfont #1: #2}%
  \ifdim \wd\@tempboxa >\hsize
    \normalsize\fontsize{10}{12}\selectfont #1: #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}

\setlength{\voffset}{0mm}
\setlength{\headheight}{25mm}
\setlength{\headsep}{0mm}
\setlength{\textheight}{250mm}
\setlength{\textwidth}{170mm}
\setlength{\textfloatsep}{12pt}
\setlength{\abovedisplayskip}{6pt}
\setlength{\belowdisplayskip}{6pt}
\setlength{\parindent}{0mm}
\setlength{\abovecaptionskip}{6pt}
\setlength{\belowcaptionskip}{6pt}
\setlength{\itemsep}{0pt}

\usepackage[paperheight=29.7cm, paperwidth=21cm,nomarginpar,textwidth=17cm ,headheight=2.5cm,top=1.5cm,bottom=1.5cm]{geometry}
\headsep=0pt

\AtBeginDocument{%
  \pagestyle{empty}
  \maketitle
  \thispagestyle{empty}
  \let\maketitle\relax
}

\endinput
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{scoreLnplot}
\alias{scoreLnplot}
\title{Generate Score Line Plot}
\usage{
scoreLnplot(
  data,
  savedir,
  pos.col = "position",
  wt.col = "wildtype",
  mean.col = "phi_mean",
  sd.col = "sigma2_mean",
  npos = 50,
  ncol = 1,
  pos.step = 1,
  x.text = 8,
  y.text = 8,
  seq.text = 3,
  ht = 30,
  wd = 15,
  name = "Lineplot",
  savepdf = TRUE,
  show = FALSE
)
}
\arguments{
\item{data}{Position data frame. Expected to have columns containing information about
position, control amino acid, mutated amino acid, mutation type, and score.}

\item{savedir}{Character string specifying the directory to save plots.}

\item{pos.col}{Column name in \code{data} for mutation positions. Default is \code{position}.}

\item{wt.col}{Column name in \code{data} for wildtype amino acids. Default is \code{wildtype}.}

\item{mean.col}{Column name in \code{data} for position-level parameters. Default is \code{phi_mean}.}

\item{sd.col}{Column name in \code{data} for error bar of mean parameters. Default is \code{sigma2_mean}.}

\item{npos}{Integer specifying the number of positions per subplot. Default is \code{50}.}

\item{ncol}{Integer specifying the number of columns of subplots. Default is \code{1}.}

\item{pos.step}{Integer specifying the steps between x-axis labels. Default is \code{5}.}

\item{x.text}{Numeric value for x-axis text size. Default is \code{8}.}

\item{y.text}{Numeric value for x-axis text size. Default is \code{8}.}

\item{seq.text}{Numeric value for wildtype sequence text size. Default is \code{3}}

\item{ht}{Numeric value for the height of the saved plot. Default is \code{30}.}

\item{wd}{Numeric value for the width of the saved plot. Default is \code{15}.}

\item{name}{Character string specifying the base name of the saved file.}

\item{savepdf}{Logical indicating whether to also save a PDF version of the plot.
Default is \code{TRUE}.}

\item{show}{Logical indicating whether or not to display the plot in the viewer.
Default is \code{FALSE}.}
}
\value{
NULL.
}
\description{
\code{scoreLn} generates a line plot for visualizing the changes of position-level
estimate across positions.
}
\examples{
\dontrun{
scoreLnplot(data = scores.pos, mean.col = "rho_mean", sd.col = "rho_sd",
            savedir = "./plots/", name = "LinePlot")
}

}

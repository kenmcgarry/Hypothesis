# parser.R
# Parser for Boolean, algebraic and arithmetic expressions; and first-order formulas; 
# and graphical plotter using igraph

library(tokenizers)
library(igraph)

## 1. Base language

## The base language permits the formation of Boolean and algebraic terms and equations.
## The basic syntax contains open and close parentheses, and no simplication is allowed.
## The atoms are the numerals, atomic letters, variables, and the truth values.
## There are eight constructors:
## : five logical ("!","&","v",">","=").
## : three algebraic ("S", "+", "*").
## "!" and "S" are unary prefix constructors.
## "&","v",">",+","*","=" are binary infix constructors.
## The constructors are single symbols.

## The syntax of the language can be described by a recursive definition, where "$v" is 
## short for applying some unary prefix constructor $ to a string v and "(u o v)" is short for 
## applying some binary infix constructor o to strings u, v.
## (Rule 1) every atom is a term.
## (Rule 2) if v is a term, then $v is a term.
## (Rule 3) if u,v are terms, then (u o v) is a term.
## (Rule 4) is u,v are terms, then (u = v) is an equation.
## (Rule 5) Nothing else is a term or an equation.
## In particular:
## Boolean expressions like "(P & (Q v !R))" count as terms (they are formulas of propostional logic).
## Algebraic expressions like "S(x + y)" count as terms.
## Algebraic expressions like "((x + Sy) = S(x + y))" count as equations (they are algebraic atomic formulas).

## The fundamental idea of the parser, given an argument string, is to 
## build a 'construction sequence' for that string by top-down parsing,
## by iterating a 'sweep' function.
## An initial processing table is constructed for the string.
## The complexity of the string in terms of occurrences of constructors is computed.
## The 'sweep' function is then applied to the table, parsing untested components
## into their further immediate components as dictacted by the production rules.
## The sweep function is iterated as many times as the complexity of the string requires.

op.paren <- "("; cl.paren <- ")"; caret <- "^";
neg <- "!"; and <- "&"; or <- "v"; cond <- ">";
succ <- "S"; plus <- "+"; prod <- "*"; eq <- "="; zero <- "0"; one <- "1";
Parentheses <- c(op.paren,cl.paren);
Constructors <- c(neg,and,or,cond,succ,plus,prod,eq);
Atoms <- c(zero,one,"TRUE","FALSE","P","Q","R","X","Y","Z","x","y","z");

## 2. Text processing: a series of functions which manipulate, cut & paste text
## The most important is 'split.infix', which splits any compound string formed
## from a binary infix constructor (e.g., "&", "+", "=") into its two parts

clean.text <- function(txt){gsub(" ", "", txt, fixed = TRUE)}
normalize.text <- function(txt){
  txt <- clean.text(txt);
  l <- length(txt);
  if (l == 0) txt <- as.character(FALSE);
  if (is.na(txt)) txt <- as.character(FALSE);
  if (txt == "") txt <- as.character(FALSE);
  if (is.null(txt)) txt <- as.character(FALSE);
  txt
}
tok <- function(txt){
  txt <- normalize.text(txt);
  tokenize_characters(txt, lowercase = FALSE,strip_non_alphanum = FALSE)[[1]]
}
rep.symbol <- function(symbol,n){
  seq <- rep(0,n+1);
  seq[1] <- "";
  if (n != 0){for (i in 1:n){seq[i+1] <- paste0(seq[i],symbol)}};
  seq[n+1]
}
extract.first.symbol <- function(txt){
  txt <- normalize.text(txt);
  X <- tok(txt);
  n <- length(X);
  Out <- FALSE
  if (txt == FALSE) Out <- NULL;
  if (txt != FALSE) Out <- X[1];
  Out
}
extract.final.symbol <- function(txt){
  txt <- normalize.text(txt);
  X <- tok(txt);
  n <- length(X);
  Out <- FALSE
  if (txt == FALSE) Out <- NULL;
  if (txt != FALSE) Out <- X[n];
  Out
}
has.open.par.prefix <- function(txt){
  txt <- extract.first.symbol(txt);
  if (is.null(txt)) FALSE else if (txt == "(") TRUE else FALSE;
}
has.close.par.suffix <- function(txt){
  txt <- extract.final.symbol(txt);
  if (is.null(txt)) FALSE else if (txt == ")") TRUE else FALSE;
}
enclose.with.parentheses <- function(txt){paste0("(",txt,")")}
is.parenthesized <- function(txt){
  P1 <- has.open.par.prefix(txt); P2 <- has.close.par.suffix(txt);
  if ((P1 == TRUE) & (P2 == TRUE))TRUE else FALSE}
deparenthesize <- function(txt){
  if (is.parenthesized(txt) == FALSE) FALSE;
  X <- tok(txt);
  n <- length(X); 
  if (n == 0)FALSE else txt.1 <- substr(txt,2,n-1);
  if (is.parenthesized(txt) == TRUE) txt.1 else FALSE}
paren.test <- function(txt){
  X <- tok(txt);
  w.open <- which(X == "(");
  w.close <- which(X == ")");
  if (length(w.open) == length(w.close))1 else 0
}
split.infix <- function(txt,split){
  k <- is.parenthesized(txt);
  X <- tok(txt);
  txt <- deparenthesize(clean.text(txt));
  n <- length(X);
  w <- which(X == split)-1;
  len <- length(w);
  if ((n == 1)|(len == 0)) Out <- "FALSE";
  T <- matrix(0,len,2);
  if (len != 0) {for (i in 1:len){T[i,1] <- substr(txt, 1,w[i]-1); T[i,2] <- substr(txt,w[i]+1,n)}};
  P <- matrix(0,len,2);
  if (len != 0 & n != 0) for (i in 1:len){P[i,1] <- paren.test(T[i,1]); P[i,2] <- paren.test(T[i,2])}
  i <- which(P[,1] == P[,2] & P[,1] == 1);
  if (length(i) == 0) Out <- FALSE else Out <- c(T[i,1],T[i,2]);
  if (k == FALSE) Out <- FALSE;
  Out
}
apply.prefix <- function(txt,prefix){clean.text(paste(prefix,txt))};
apply.suffix <- function(txt,suffix){clean.text(paste(txt,suffix))};
has.prefix <- function(txt,prefix){
  X <- tok(txt);
  n <- length(X); 
  if (n > 0) f <- X[1] else f <- 0;
  if (f == prefix){TRUE} else FALSE
}
delete.prefix <- function(txt,prefix){
  X <- tok(txt);
  n <- length(X); 
  if (has.prefix(txt,prefix) == TRUE) {substr(clean.text(txt),2,n)} else FALSE
}
has.suffix <- function(txt,suffix){
  X <- tok(txt);
  n <- length(X); 
  if (n > 0) f <- X[n] else f <- 0;
  if (f == suffix){TRUE} else FALSE
}
delete.suffix <- function(txt,suffix){
  X <- tok(txt);
  n <- length(X); 
  if (has.suffix(txt,suffix) == TRUE) {substr(clean.text(txt),n-1,n)} else FALSE
}
has.unary.prefix <- function(txt,symbol){
  X <- tok(txt);
  n <- length(X); 
  if (n > 0) f <- X[1] else f <- 0;
  if (f == symbol){TRUE} else FALSE
}
delete.unary.prefix <- function(txt,symbol){
  X <- tok(txt);
  n <- length(X); 
  if (has.unary.prefix(txt,symbol) == TRUE) {substr(clean.text(txt),2,n)} else FALSE
}
apply.infix <- function(txt1,txt2,infix){
  clean.text(enclose.with.parentheses(paste(txt1,C,txt2)))
}
prefix.numeral <- function(x,root,prefix){
  if(x==0)
    return(root)
  else
    return(apply.prefix(prefix.numeral(x-1,root,prefix),prefix))
}
suffix.numeral <- function(x,root,suffix){
  if(x==0)
    return(root)
  else
    return(apply.suffix(suffix.numeral(x-1,root,suffix),suffix))
}
count.occurrences <- function(txt,symbol){t <- tok(txt); w <- which(t == symbol); length(w)}

## 3. Atoms, literals, constructors, variables, numerals.

## Numerals are formed by prefixing n (>=0) occurrences of "S" to "0" (e.g., "SSSS0").
## Extra variables are formed by suffixing n (>=0)occurrences of "^" to "x" (e.g., "x^^").
## There are a series of tests and functions which check whether the string
## is a primitive, an atom, a constructor, a literal, a numeral, etc.
## The syntactic operators corresponding to the constructors are also defined.

is.primitive <- function(txt,X){any(X == txt)}
is.atom <- function(txt){is.primitive(txt,Atoms)|is.variable(txt)}
is.neg.atom <- function(txt){
  has.unary.prefix(txt,"!") & is.atom(delete.unary.prefix(txt,"!"))
}
is.literal <- function(txt){
  is.atom(txt)|is.neg.atom(txt)
}
numeral <- function(x){
  if(x==0)
    return("0")
  else
    return(SUCC(numeral(x-1)))
}
is.numeral <- function(txt){
  txt <- normalize.text(txt);
  X <- tok(txt);
  n <- length(X)+2
  T <- rep(0,n)
  for (i in 1:n){T[i] <- numeral(i)}
  any(T == txt)
}
variable <- function(x){
  if(x==0)
    return("x")
  else
    return(clean.text(paste(variable(x-1),"^")))
}

is.variable <- function(txt){
  txt <- normalize.text(txt);
  X <- tok(txt);
  n <- length(X)+2
  T <- rep(0,n)
  for (i in 1:n){T[i] <- variable(i)}
  any(T == txt)
}

is.constructor <- function(txt){is.primitive(txt,Constructors)};
apply.unary.constructor <- function(txt,C){clean.text(paste(C,txt))};
apply.binary.constructor <- function(txt1,txt2,C){
  clean.text(enclose.with.parentheses(paste(txt1,C,txt2)))
}

NEG <- function(txt){apply.unary.constructor(txt,"!")}
SUCC <- function(txt){apply.unary.constructor(txt,"S")}
AND <- function(txt1,txt2){apply.binary.constructor(txt1,txt2,"&")}
OR <- function(txt1,txt2){apply.binary.constructor(txt1,txt2,"v")}
IF <- function(txt1,txt2){apply.binary.constructor(txt1,txt2,">")}
PLUS <- function(txt1,txt2){apply.binary.constructor(txt1,txt2,"+")}
PROD <- function(txt1,txt2){apply.binary.constructor(txt1,txt2,"*")}
EQ <- function(txt1,txt2){apply.binary.constructor(txt1,txt2,"=")}

## 4. Compound expressions: complexity, tests and processing matrices

## The first six functions compute complexity of the string in terms
## presence of logical and algbraic constructors
## The third is a test which the parser used to see if the string is possibly
## a binary compund
## The final two are processing matrices used to build up the construction sequence;
## They take a string and a constructor, and output its immediate components,
## irrespective of whether they're well-formed

count.occurrences <- function(txt,symbol){t <- tok(txt); w <- which(t == symbol); length(w)}
unary.log.complexity <- function(txt){count.occurrences(txt,"!")}
unary.alg.complexity <- function(txt){count.occurrences(txt,"S")}
binary.log.complexity <- function(txt){count.occurrences(txt,"v") + count.occurrences(txt,"&") + count.occurrences(txt,">")}
binary.alg.complexity <- function(txt){count.occurrences(txt,"+") + count.occurrences(txt,"prod") + count.occurrences(txt,"=")}
complexity <- function(txt){unary.log.complexity(txt) + unary.alg.complexity(txt) + binary.log.complexity(txt) + binary.alg.complexity(txt)}
test.binary.symb <- function(txt,symb){count.occurrences(txt,symb)>0 & is.parenthesized(txt) & split.infix(txt,symb)[1] != FALSE}
unary.matrix <- function(txt,C){rbind(delete.unary.prefix(txt,C),0,txt)}
binary.matrix <- function(txt,C){rbind(c(split.infix(txt,C)[2],split.infix(txt,C)[1]),c(0,0),c(txt,txt))}

## 5. Parser algorithm

## The main engine of the parser is the 'Sweep' function which processes and extends
## a given table D keeping trakc of which strings have been processed.
## The parser takes a string, measures its complexity n, and bulds the initial processing table.
## It then iterates the 'sweep' function n times. Its output is the finished table.
## Reducer is a function which checks that the construction sequence terminates with atoms
## and that the original string is built from the constructors as required by the
## production rules.
## The 'Construction.sequence' function produces a table of the completed sequence.

Sweep <- function(D){
  w <- which(D[2,] == 0); 
  Component <- function(i){D[1,w[i]]};
  len <- length(w);
  C <- rep(0,len);
  for (i in 1:len){C[i] <- Component(i)};
  for (i in 1:len){
    t <- C[i];
    D[2,w[i]] <- 1;
    if (has.unary.prefix(t,"!") == TRUE) {D <- cbind(D,rbind(unary.matrix(t,"!"),w[i],"!"))};
    if (has.unary.prefix(t,"S") == TRUE) {D <- cbind(D,rbind(unary.matrix(t,"S"),w[i],"S"))};
    if (test.binary.symb(t,"v") == TRUE) {D <- cbind(D,rbind(binary.matrix(t,"v"),w[i],"v"))};
    if (test.binary.symb(t,"&") == TRUE) {D <- cbind(D,rbind(binary.matrix(t,"&"),w[i],"&"))};
    if (test.binary.symb(t,">") == TRUE) {D <- cbind(D,rbind(binary.matrix(t,">"),w[i],">"))};
    if (test.binary.symb(t,"+") == TRUE) {D <- cbind(D,rbind(binary.matrix(t,"+"),w[i],"+"))};
    if (test.binary.symb(t,"*") == TRUE) {D <- cbind(D,rbind(binary.matrix(t,"*"),w[i],"*"))};
    if (test.binary.symb(t,"=") == TRUE) {D <- cbind(D,rbind(binary.matrix(t,"="),w[i],"="))};
  };
  D
}

Parser <- function(txt){
  txt <- clean.text(txt);
  Vector <- c(txt);
  Status <- c(0);
  Predecessors <- c(0);
  Pred.label <- c(0);
  Connect <- c(0);
  Table <- rbind(Vector, Status, Predecessors, Pred.label, Connect);
  Comp <- complexity(txt);
  for (j in 1:(Comp+1)){Table <- Sweep(Table)};
  X <- as.vector(rev(Table[1,]));
  N <- length(X);
  Y <- as.vector(rev(Table[3,]));
  W <- as.vector(rev(Table[5,]));
  Z <- rep(0,N)
  if (N > 1){Z <- (N+1) - as.numeric(rev(Table[4,]))};
  if (N == 1){Z <- c(1)};
  Y[N] <- "Root";
  Z[N] <- "NULL";
  M <- as.data.frame(cbind(X,Y,Z,W));
  colnames(M) <- c("Construction", "Parent", "Parent Index","Constructor");
  M
}

Reducer <- function(C,A){
  N <- length(C);
  Reduced <- rep(0,N); 
  Reduced <- A;
  for (i in 1:N) {
    M <- as.vector(Parser(C[i])[,3]);
    if (M[1] == "NULL") Reduced[i] <- FALSE;
    if (M[1] != "NULL"|A[i]==TRUE) Reduced[i] <- TRUE}
  Reduced
}

Construction.sequence <- function(txt){
  C <- as.vector(Parser(txt)[,1]);
  index <- as.vector(Parser(txt)[,3]);
  conn <- as.vector(Parser(txt)[,4]);
  N <- length(C);
  AT <- rep(0,N); T <- rep(0,N); conn[N] <- "NULL";
  for (i in 1:N){if (is.atom(C[i])) AT[i] <- TRUE};
  B <- as.data.frame(cbind(C,AT,index,conn));
  colnames(B) <- c("Construction", "Atom","Parent Index","Constructor");
  B
}

## 6. Syntactical analysis

## These functions perform tests on the input string, identifying if it is wellformed, 
## what atoms occur, ## what the main constructor is, what the immediate precessors are.
## The fifth function 'Syntax.test' summarizes some of this data.

is.wellformed <- function(txt){
  C <- as.vector(Parser(txt)[,1]);
  P <- as.vector(Parser(txt)[,2]);
  L <- as.vector(Parser(txt)[,3]);
  N <- length(C);
  AT <- rep(0,N); T <- rep(0,N); P[N] <- "Root";
  for (i in 1:N){if (is.atom(C[i])) AT[i] <- TRUE};
  Reduced <- Reducer(C,AT);
  if (any(Reduced == FALSE)) Out <- FALSE else Out <- TRUE;
  Out
}

atoms.in.expression <- function(txt){
  T <- Construction.sequence(txt);
  I <- T[,2];
  w <- which(I == 1);
  X <- as.character(T[,1]);
  n <- length(w);
  a <- rep(0,n);
  if (n != 0) {for (i in 1:n) {a[i] <- as.character(X[w[i]])}}
  if (n != 0) Out <- a else Out <- NULL;
  rev(Out)
}

main.constructor <- function(txt){
  f <- is.wellformed(txt)
  a <- is.atom(txt)
  T <- Construction.sequence(txt)
  C <- as.character(T[,4])
  n <- length(C)
  if (f == TRUE & a == TRUE){Out <- "Atom"}
  if (f == TRUE & a == FALSE){Out <- C[n-1]}
  if (f == FALSE){Out <- NULL}
  Out
}

imm.pred <- function(txt){
  T <- as.matrix(Construction.sequence(txt))
  forms <- as.vector(T[,1]);
  indices <- as.vector(T[,3]);
  a <- is.wellformed(txt);
  n <- length(indices);
  w <- which(indices == n);
  F <- forms[w];
  k <- length(w);
  Out <- FALSE
  if (a == TRUE & k == 1) Out <- F[1];
  if (a == TRUE & k == 2) Out <- rev(F[1:2]);
  rev(Out)
}

is.term <- function(txt){is.wellformed(txt) & main.constructor(txt) != "="}

is.equation <- function(txt){is.wellformed(txt) & main.constructor(txt) == "="}

syntax.test <- function(txt){
  Table <- matrix(0,3,3);
  F <- is.wellformed(txt);
  M <- main.constructor(txt)
  I <- imm.pred(txt);
  num.pred <- length(I)
  A <- atoms.in.expression(txt);
  num.atoms <- length(A)
  p1 <- 4 + num.pred - 1;
  p2 <- p1 + 1;
  p3 <- p2 + num.atoms -1;
  Output <- c(c(clean.text(txt)),F,M,I,A);
  n <- length(Output);
  Data <- rep(0,n)
  Data[1] <- "Input Expression";
  Data[2] <- "Is Wellformed?";
  Data[3] <- "Main Constructor";
  for (i in 4:p1){Data[i] <- "Immediate Predecessor"}
  for (i in p2:p3){Data[i] <- "Atom"}
  if (F == FALSE) Data <- Data[1:3];
  if (F == FALSE) Output <- Output[1:3];
  as.data.frame(cbind(Data,Output))
}

## 7. Parsing tree graphical plotters

## These two functions plot a parsing tree for a string.
## The first takes a string as argument;d the second takes a construction sequence as argument.

Plot.parse.tree <- function(string){
  CS <- Construction.sequence(string)
  text <- as.character(CS$Construction);
  N <- length(text);
  parent <- as.vector(CS$`Parent Index`);
  Nodes <- 1:N
  to <- Nodes[1:(N-1)]
  from <- as.numeric(parent[1:(N-1)])
  Rel <- as.data.frame(cbind(from,to))
  g <- graph_from_data_frame(Rel, directed=TRUE, vertices=Nodes)
  labels <- as.character(CS$Constructor)[1:(N-1)]
  V(g)$label <- text
  E(g)$label <- labels
  plot(g, layout=layout.reingold.tilford,vertex.shape = "circle", 
       vertex.size=2, vertex.color = "yellow", 
       vertex.label.degree=0, vertex.label.dist=1.2,
       edge.arrow.size=0.1, edge.color = "green",
       main = "Parsing Tree")
}
Plot.parse.tree.cs <- function(CS){
  text <- as.character(CS$Construction);
  N <- length(text);
  parent <- as.vector(CS$`Parent Index`);
  Nodes <- 1:N
  to <- Nodes[1:(N-1)]
  from <- as.numeric(parent[1:(N-1)])
  Rel <- as.data.frame(cbind(from,to))
  g <- graph_from_data_frame(Rel, directed=TRUE, vertices=Nodes)
  labels <- as.character(CS$Constructor)[1:(N-1)]
  V(g)$label <- text
  E(g)$label <- labels
  plot(g, layout=layout.reingold.tilford,vertex.shape = "circle", 
       vertex.size=2, vertex.color = "yellow", 
       vertex.label.degree=0, vertex.label.dist=1.2,
       edge.arrow.size=0.1, edge.color = "green",
       main = "Parsing Tree")
}

## 8. Some examples

Parser(PLUS(numeral(8),zero))
Construction.sequence("!!(P v (P & (TRUE v FALSE)))")
Construction.sequence("!((x^^ + y) = SSSSSSSSSSx)")
Construction.sequence("((T + F) = (F + T))")
Construction.sequence(PLUS(SUCC(variable(0)),variable(4)))
Construction.sequence(PLUS(numeral(8),zero))
is.wellformed("P & Q")
is.wellformed("(P & Q")
is.wellformed("(P & Q)")
is.term("(x^^ + y)")
is.equation("((x^^ + y) = x)")
main.constructor("!!P")
main.constructor("((x + y) = (y + x))")
imm.pred("((x + y) = (y + x))")
syntax.test("((x^^ + y) = x)")
syntax.test("SSx^^")
syntax.test("((x + y) = SSSSSSSSSSx)")
syntax.test("")
Plot.parse.tree("S(x * SSx^^)")
Plot.parse.tree("!!(P > (P v P))")
Plot.parse.tree("!!(P v (Q & (TRUE v FALSE)))")
T <- Construction.sequence(PLUS(numeral(8),zero)); 
Plot.parse.tree.cs(T)
string <- EQ(PROD("x",SUCC("y")),PLUS(PROD("x","y"),"x")); 
Plot.parse.tree(string)
T <- Construction.sequence(PLUS(numeral(8),zero));
Plot.parse.tree.cs(T)
string <- "!((Sx = 0) > (P > (Q > !!!P)))"
Construction.sequence(string)
Plot.parse.tree(string)




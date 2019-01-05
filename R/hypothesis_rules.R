# hypothesis_rules.R
# Hypothesis (theory) revision using Bayesian techniques.

# https://projecteuclid.org/download/pdf_1/euclid.ndjfl/1093636527
# http://www.logicinaction.org/docs/ch3.pdf
# http://www.logicinaction.org/docs/ch4.pdf
# https://www.dcode.fr/syllogism

# There are an infinity of syllogisms, but structurally there are only 24 types of 
# syllogisms which can have a logical conclusion. The conclusion is a purely logical 
# conclusion. If at least one of the 2 assertions / premises is false, then it is possible 
# to deduce things that are completely false, although perfectly logical. Similarly it is 
# possible to deduce something right from false or incomplete premises.
# For a syllogism to be valid, it must consist of two premises. A premise is itself 
# made up of a quantity (All, None, etc.), a subject, a verb (is / is not) and a predicate 
# (a subject attribute). The two premises must also have a middle subject or predicate, ie it 
# must be part of the 2 premises.
# Example: No A are B. All C are A. => No C are B.

# In order to be a standard-form categorical syllogism, three requirements must be met:
#  1) All three statements must be standard-form categorical propositions.
#  2) The two occurrences of each term must be identical and have the same sense.
#  3) The major premise must occur first, the minor premise second, and the conclusion last.
# The mood of a categorical syllogism consists of the type of categorical propositions 
# involved (A, E, I, or O) and the order in which they occur. The middle term can be arranged in 
# the two premises in four different ways. These placements determine the figure of the categorical syllogism.
# There are six rules for standard-form categorical syllogisms:
#  1) The middle term must be distributed in at least one premise.
#  2) If a term is distributed in the conclusion, then it must be distributed in a premise.
#  3) A categorical syllogism cannot have two negative premises.
#  4) A negative premise must have a negative conclusion.
#  5) A negative conclusion must have a negative premise.
#  6) Two universal premises cannot have a particular conclusion.


# https://www.rdocumentation.org/packages/NLP/versions/0.2-0/topics/Tree
p <- Tree("VP", list(Tree("V",list("saw")),
               Tree("NP",
                    list("him"))))

p <- Tree("S", list(Tree("NP",
                    list("I")),p))
p

print(p, width = 10)

s <- "(S (NP I) (VP (V saw) (NP him)))"  # this is how users will create their syllogistic statements
p1 <- Tree_parse(s)  # Tree_parse() converts the bracket string into tree structure
print(p1,width=10)

## Extract the leaves by recursively traversing the children and
## recording the non-tree ones:
Tree_leaf_gatherer <-
  function()
  {
    v <- list()
    list(update =
           function(e) if(!inherits(e, "Tree")) v <<- c(v, list(e)),
         value = function() v,
         reset = function() { v <<- list() })
  }
g <- Tree_leaf_gatherer()
y <- Tree_apply(p, g$update, recursive = TRUE)
g$value()


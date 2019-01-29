library("mstate")
data("ebmt4")
ebmt <- ebmt4
head(ebmt)


tmat <- transMat(x = list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6), c(5, 6),c(), c()), names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death"))
tmat


msebmt <- msprep(data = ebmt, trans = tmat, time = c(NA, "rec", "ae", "recae", "rel", "srv"), status = c(NA, "rec.s", "ae.s", "recae.s", "rel.s", "srv.s"), keep = c("match", "proph", "year", "agecl"))

macro_rules! i {
    ($v:ident[$($x:expr),*]) => {
        &$v.slice(s![$($x),*])
    };
}

pub(crate) use i;

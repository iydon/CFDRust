macro_rules! i {
    ($a:ident[$($x:expr),*]) => {
        &$a.slice(s![$($x),*])
    };
}

macro_rules! assign {
    ($a:ident[$($x:expr),*] = $rhs:expr) => {
        $a.slice_mut(s![$($x),*]).assign(&$rhs)
    }
}

macro_rules! fill {
    ($a:ident[$($x:expr),*] = $v:expr) => {
        $a.slice_mut(s![$($x),*]).fill($v)
    }
}

pub(crate) use {assign, fill, i};

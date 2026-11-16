ACCTS=(<account> <account> <account> <account>)
EXCL='gnode001,gnode002,gnode003,gnode004,gnode005,gnode006,gnode007,gnode008,gnode011,gnode012,gnode017,gnode020,gnode021,gnode022,gnode023,gnode025,gnode026,gnode027,gnode028,gnode029,gnode031,gnode033,gnode034,gnode035,gnode036,gnode038,gnode040,gnode041,gnode042,gnode043,gnode044,gnode066,gnode077'
MASK=':392-406'
declare -A AJOB=( [pb_bhd2]=48756 [docked_bhd2]=48757 )
declare -A STATE=( [pb_bhd2]=pb [docked_bhd2]=docked )
SSH="ssh -o BatchMode=yes -o ConnectTimeout=20 -o ControlMaster=auto -o ControlPath=/tmp/cm-%r -o ControlPersist=10m"
# submit window $3 of leg $1 as account $2, optional dependency $4
submit_win() {
  local leg=$1 acct=$2 win=$3 dep=${4:-}
}
